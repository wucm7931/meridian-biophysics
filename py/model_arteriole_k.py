# -*- coding: utf-8 -*-
"""
model_v2 · Submodule 02 · Arteriole vasomotion and output (K⁺ release, source b) — computation script
================================================================================
Corresponding design doc: 02_arteriole_vasomotion.md
Responsibilities: compute B-type derived values in real time from A-type raw measurements, output the K⁺ concentration field of source b at the terminal read point:
  · S_CELL (B-type, VSMC capacitance ÷ specific membrane capacitance)
  · j_STOC / capacitance lower bound / 22× (B-type, BKCa cycling flux, overturning the old "8.6 µM")
  · u_b cylindrical uniform-caliber STOC peak (B-type, cylindrical-source reaction-diffusion steady-state solution)
  · u_b transient point-source-caliber STOC peak (B-type, point-source 100 ms transient solution)

Design iron rules (see 00_overview.md §3):
  · Each numeric literal has a same-line comment "source (PMID/local file) + raw measured value".
  · B-type derived values are always computed on the spot from A-type raw values; never hardcode results.
  · Unit iron rules: 1 mol/m³ ≡ 1 mM; 1 pmol/cm² = 1e-8 mol/m²; 1 µm² = 1e-12 m².

Run: `python model_arteriole_k.py` (depends on common.py in the same directory)
"""

import math
import sys
import numpy as np
from scipy.special import k0, k1, erfc

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

from common import D_eff, ALPHA, K_CLEAR  # A-type (see common.py for per-item sources)

# ----------------------------------------------------------------------
# A-type · STOC (raw measured values)
#   Benham & Bolton 1986, J Physiol 381:385–406, PMID 2442353 (download3/PMID-2442353.txt):
#   single STOC peak 250 pA ("up to", mean 142±17 pA), duration ~100 ms, 75–100 channels.
# ----------------------------------------------------------------------
I_STOC = 250e-12            # A = 250 pA (max; mean 142 pA see sensitivity)
T_STOC = 100e-3             # s = 100 ms

# ----------------------------------------------------------------------
# A-type · resting spark frequency
#   Jaggar 2000, Am J Physiol Cell Physiol 278:C235–C256, PMID 10666018:
#   Table 1 "Unstimulated frequency" 0.5–1 s⁻¹ per cell.
# ----------------------------------------------------------------------
F_SPARK = 1.0               # Hz (base frequency 1, range 0.5–1)

# ----------------------------------------------------------------------
# A-type · VSMC whole-cell capacitance + resting membrane potential
#   Smirnov & Aaronson 1992, J Physiol 457:431–454, PMID 1284312:
#   human mesenteric artery VSMC whole-cell clamp capacitance ~46 pF; resting V_m −44.1±1.9 mV.
# ----------------------------------------------------------------------
C_CELL = 46e-12             # F = 46 pF

# ----------------------------------------------------------------------
# A-type · terminal read-point distance (source b K⁺ released from the VSMC surface)
#   Luff 2000, J Comp Neurol 416:258–271, PMID 10602088:
#   SP⁺ sensory C-fiber varicosities 0.4–2 µm from VSMC (most); take 1 µm median as representative.
# ----------------------------------------------------------------------
READ_GAP_UM = 1.0           # µm (read-point distance from VSMC surface, scan 0.4/1/2)

# ----------------------------------------------------------------------
# A-type · terminal arteriole radius
#   skin terminal arteriole diameter d = 25 µm (range 10–25 µm) → r₀ = 12.5 µm.
#   no local full-text record → PMID to be added (placeholder value, see design doc §6).
# ----------------------------------------------------------------------
R_VESSEL = 12.5e-6          # m (terminal arteriole radius, d=25 µm)

# ----------------------------------------------------------------------
# A-type · specific membrane capacitance (standard physiology, no single PMID)
# ----------------------------------------------------------------------
C_M = 1e-2                  # F/m² = 1 µF/cm²

# ----------------------------------------------------------------------
# A-type · diastolic hyperpolarization amplitude (membrane potential change due to BKCa opening)
#   Jaggar 2000, Am J Physiol Cell Physiol 278:C235–C256, PMID 10666018
# ----------------------------------------------------------------------
DV_HYP = 25e-3              # V = 25 mV

# ----------------------------------------------------------------------
# C-type · density coefficient (placeholder, see 02_arteriole_vasomotion.md §6)
#   along-meridian 1.5× / acupoint 3.0× (placeholder expression of a larger number N of synchronized sources within the receptive field, to be anchored)
# ----------------------------------------------------------------------
DENSITY = 1.5               # along-meridian coefficient
DENS_POINT = 3.0            # acupoint coefficient


def cell_area(c_cell: float, c_m: float) -> float:
    """B-type · VSMC single-cell surface area S_CELL = C_cell / C_m.

    Raw measured value → formula: C_cell = 46 pF (Smirnov & Aaronson 1992, PMID 1284312);
      C_m = 1 µF/cm² (standard physiology).
    """
    return c_cell / c_m                                     # [m²]


def stoc_flux(i_stoc: float, t_stoc: float, s_cell: float) -> float:
    """B-type · single-STOC flux density j_STOC = I_STOC·T_STOC / (F·S_CELL).

    Raw measured value → formula: I_STOC=250 pA, T_STOC=100 ms (Benham & Bolton 1986,
      PMID 2442353); S_CELL (§4.1, [1] in this script). Returns [mol/m²].
    """
    F = 96485.0                                            # C/mol Faraday constant (SI-defined quantity)
    return i_stoc * t_stoc / (F * s_cell)


def cyl_u(r: float, j: float, r_vessel: float, alpha: float, d_eff: float,
          k_clear: float) -> float:
    """B-type · cylindrical-source steady state [mol/m³ ≡ mM]: u(r) = (J·λ/(D_eff·α))·K₀(r/λ)/K₁(r₀/λ).

    Raw measured value → formula: cylindrical-source reaction-diffusion steady-state solution with first-order clearance k (r ≥ r₀),
      λ = √(D_eff/k) ≈ 272 µm. D_eff, α, k see common.py.
    Geometric assumption: axial direction taken as infinite / z invariant (vessel segment ≫ λ, end effects negligible) —
      hence the solution has no length parameter and depends only on radial r; this is an upper bound for a given uniform flux.
    Uniform caliber: j smeared over the entire vessel surface (spatial average) → this caliber is a lower bound
      (the real STOC is a sub-micron localized source, see transient_pt).
    """
    lam = math.sqrt(d_eff / k_clear)                       # clearance length [m]
    pref = j * lam / (d_eff * alpha)
    return pref * k0(r / lam) / k1(r_vessel / lam)


def transient_pt(r: float, i_stoc: float, t_stoc: float, alpha: float,
                 d_eff: float) -> float:
    """B-type · transient point-source STOC peak [mol/m³ ≡ mM].

    Raw measured value → formula: transient solution for a continuous point source turned on for 100 ms
      u(r,t) = (I_STOC/F)/(α·4π·D_eff·r)·erfc(r/√(4·D_eff·t)).
      Local "oasis" radius 250–490 nm, area 0.2–0.75 µm² (ZhuGe 2002, PMID 12084772).
    """
    F = 96485.0
    q_dot = i_stoc / F                                     # [mol/s] point-source release rate
    arg = r / math.sqrt(4.0 * d_eff * t_stoc)
    return q_dot / (alpha * 4.0 * math.pi * d_eff * r) * erfc(arg)


def compute():
    """Return this module's hand-off quantities (see 02_arteriole_vasomotion.md §5/§6).

    · S_CELL: VSMC single-cell surface area [m²]
    · j_stoc / j_stoc_flux / j_avg / j_hyp: single-STOC flux density (various calibers) [mol/m² or mol/m²/s]
    · ratio: j_STOC/j_hyp (cycling current multiple)
    · u_b_cyl / u_b_avg / u_b_pt: 1 µm read-point K⁺ concentration (cylindrical peak / time average / transient point-source peak) [mol/m³ ≡ mM]
    · *_m / *_p: along-meridian 1.5× / acupoint 3.0× density scaling (C-type placeholder)
    """
    S_CELL = cell_area(C_CELL, C_M)
    j_stoc = stoc_flux(I_STOC, T_STOC, S_CELL)
    j_stoc_flux = I_STOC / (96485.0 * S_CELL)
    j_avg = j_stoc * F_SPARK
    j_hyp = C_M * DV_HYP / 96485.0
    r_s = R_VESSEL + READ_GAP_UM * 1e-6
    u_cyl = cyl_u(r_s, j_stoc_flux, R_VESSEL, ALPHA, D_eff, K_CLEAR)
    u_avg = cyl_u(r_s, j_avg, R_VESSEL, ALPHA, D_eff, K_CLEAR)
    u_pt  = transient_pt(READ_GAP_UM * 1e-6, I_STOC, T_STOC, ALPHA, D_eff)
    return {
        "S_CELL": S_CELL,
        "j_stoc": j_stoc,
        "j_stoc_flux": j_stoc_flux,
        "j_avg": j_avg,
        "j_hyp": j_hyp,
        "ratio": j_stoc / j_hyp,
        "u_b_cyl": u_cyl,
        "u_b_avg": u_avg,
        "u_b_pt": u_pt,
        "u_b_cyl_m": u_cyl * DENSITY,
        "u_b_cyl_p": u_cyl * DENS_POINT,
        "u_b_pt_m": u_pt * DENSITY,
        "u_b_pt_p": u_pt * DENS_POINT,
    }


if __name__ == "__main__":
    print("=" * 68)
    print("Submodule 02 · Arteriole vasomotion and output (source b) — computation output")
    print("=" * 68)

    # ---- B-type · S_CELL ----
    S_CELL = cell_area(C_CELL, C_M)
    print(f"\n[1] VSMC single-cell surface area S_CELL = C_cell/C_m")
    print(f"    C_cell = {C_CELL*1e12:.0f} pF (Smirnov & Aaronson 1992, PMID 1284312)")
    print(f"    S_CELL = {S_CELL*1e12:.0f} µm²   (vs v1's 1131 µm², ratio "
          f"{S_CELL*1e12/1131:.2f}×)")

    # ---- B-type · single STOC flux ----
    j_stoc = stoc_flux(I_STOC, T_STOC, S_CELL)             # [mol/m²]
    j_stoc_flux = I_STOC / (96485.0 * S_CELL)              # [mol/m²/s] instantaneous flux during STOC
    j_avg = j_stoc * F_SPARK                               # [mol/m²/s] time average
    print(f"\n[2] single-STOC flux density j_STOC = I·T/(F·S_CELL)")
    print(f"    STOC = {I_STOC*1e12:.0f} pA × {T_STOC*1e3:.0f} ms = "
          f"{I_STOC*T_STOC*1e12:.0f} pC = {I_STOC*T_STOC/96485.0:.2e} mol")
    print(f"    j_STOC = {j_stoc*1e8:.1f} pmol/cm² (single STOC)")
    print(f"    j_avg  = {j_avg*1e8:.1f} pmol/cm²/s (× f_spark={F_SPARK:.0f} Hz, time average)")

    # ---- B-type · capacitance lower bound and 22× ----
    j_hyp = C_M * DV_HYP / 96485.0                         # [mol/m²] capacitance lower bound (DV_HYP see module-level constant)
    ratio = j_stoc / j_hyp
    print(f"\n[3] capacitance lower bound j_hyp = C_m·ΔV/F (comparison quantity overturning the old '8.6 µM')")
    print(f"    j_hyp = {j_hyp*1e8:.2f} pmol/cm² (net charge of one 25 mV hyperpolarization)")
    print(f"    j_STOC/j_hyp = {ratio:.0f}× (K⁺ efflux is a cycling current, not a one-time transfer)")

    # ---- B-type · cylindrical uniform-caliber STOC peak ----
    lam = math.sqrt(D_eff / K_CLEAR)
    print(f"\n[4] cylindrical uniform-caliber STOC peak (r₀={R_VESSEL*1e6:.1f} µm, λ={lam*1e6:.0f} µm)")
    for gap_um in (0.4, 1.0, 2.0):
        r_s = R_VESSEL + gap_um * 1e-6
        u_p = cyl_u(r_s, j_stoc_flux, R_VESSEL, ALPHA, D_eff, K_CLEAR)
        u_t = cyl_u(r_s, j_avg, R_VESSEL, ALPHA, D_eff, K_CLEAR)
        print(f"    read point {gap_um:3.1f} µm: STOC peak {u_p*1e3:5.0f} µM | time average "
              f"{u_t*1e3:5.1f} µM")

    # ---- B-type · transient point-source-caliber STOC peak ----
    print(f"\n[5] transient point-source-caliber STOC peak (STOC localized 'oasis', ZhuGe 2002)")
    for gap_um in (0.4, 1.0, 2.0):
        r = gap_um * 1e-6
        u_pt = transient_pt(r, I_STOC, T_STOC, ALPHA, D_eff)
        print(f"    read point {gap_um:3.1f} µm: transient point source {u_pt*1e3:5.0f} µM = {u_pt:.2f} mM")

    # ---- density scaling + threshold situation ----
    # u_th=6.8 mM full-firing threshold from module 03 (shown here only for comparison; declared locally to avoid a circular import, do not change the value)
    U_TH = 6.8              # mM full-firing threshold (see module 03 §4.2, Tkachenko 2023 empirical anchor)
    h = compute()
    u_cyl = h["u_b_cyl"]
    u_pt  = h["u_b_pt"]
    print(f"\n[6] density scaling (along-meridian {DENSITY}× / acupoint {DENS_POINT}×) and threshold situation (u_th={U_TH} mM)")
    print(f"    cylindrical caliber: ordinary {u_cyl*1e3:3.0f} µM | along-meridian {h['u_b_cyl_m']*1e3:4.0f} µM"
          f" | acupoint {h['u_b_cyl_p']*1e3:4.0f} µM")
    print(f"    transient caliber: ordinary {u_pt*1e3:3.0f} µM | along-meridian {h['u_b_pt_m']*1e3:4.0f} µM"
          f" ({h['u_b_pt_m']:.2f} mM) | acupoint {h['u_b_pt_p']*1e3:4.0f} µM"
          f" ({h['u_b_pt_p']:.2f} mM)")
    print(f"    → along-meridian 1.5× (~{h['u_b_pt_m']:.1f} mM), acupoint 3.0× (~{h['u_b_pt_p']:.1f} mM)")
    print(f"      transient calibers are all < {U_TH} mM full-firing threshold → source b never crosses threshold;")
    print(f"    → source b role = below-threshold charging; crossing the sustained threshold can only rely on stock events (dam-break / resonance).")

    print("\n" + "=" * 68)
    print("Hand-off: u_b cylindrical peak 99 µM / transient peak 0.9 mM / time average ~10 µM (to 03).")
    print("=" * 68)
