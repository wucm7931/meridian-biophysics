# -*- coding: utf-8 -*-
"""
model_v2 · Submodule 01 · Sympathetic recruitment and output (computation script)
====================================================
Corresponding design doc: 01_sympathetic_recruitment.md
Responsibilities: compute B-type derived values in real time from A-type raw measurements, output this module's hand-off quantities:
  · single varicosity surface area S_v (B-type, from varicosity diameter/length → prolate spheroid formula)
  · single varicosity per-AP K⁺ efflux Q_b (B-type, efflux × S_v)
  · source a near-field instantaneous peak u_a_peak (B-type, cleft quasi-steady-state solution)

Design iron rules (see 00_overview.md §3):
  · Each numeric literal has a same-line comment "source (PMID/local file) + raw measured value".
  · B-type derived values are always computed on the spot from A-type raw values; never hardcode results.
  · Unit iron rules: 1 mol/m³ ≡ 1 mM; 1 pmol/cm² = 1e-8 mol/m²; 1 µm² = 1e-12 m².

Run: `python model_sympathetic.py` (depends on common.py in the same directory)
"""

import math
import sys

# Windows console defaults to GBK, which cannot encode characters like µ/²/×; force stdout to use UTF-8.
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

from common import D_free, N_A  # A-type · D_free free-water KCl diffusion coefficient (Syková & Nicholson 2008); N_A Avogadro constant (SI-defined)

# ----------------------------------------------------------------------
# A-type · varicosity geometry (raw measured values)
#   Goyal & Chaudhury 2013, Auton Neurosci 176:11–31, PMID 23535140:
#   sympathetic varicosities are spindle-shaped, diameter d_v = 0.5–2.0 µm, length l_v = 1–3 µm.
# ----------------------------------------------------------------------
D_V_MIN, D_V_MED, D_V_MAX = 0.5e-6, 1.0e-6, 2.0e-6   # diameter [m], median 1.0 µm
L_V_MIN, L_V_MED, L_V_MAX = 1.0e-6, 2.0e-6, 3.0e-6   # length [m], median 2.0 µm

# ----------------------------------------------------------------------
# A⚠️-type · axonal per-AP K⁺ efflux (raw measured value + trunk→varicosity extrapolation assumption)
#   Keynes & Ritchie 1965, J Physiol 179:333–367, PMID 5853895:
#   mammalian unmyelinated nerve trunk (rabbit vagus C fibers) per-impulse K⁺ efflux = 1.0 pmol/cm².
#   ⚠️ Assumption: this value was measured on the nerve trunk and is carried over here to the varicosity membrane — same axon, same membrane, same AP;
#      identical per-area efflux is a first-order assumption, not independently verified on varicosity swelling.
#   Upper-bound check: charge conservation C_m·ΔV/F = 1.04 pmol/cm² (see [3]), independent of channel type.
# ----------------------------------------------------------------------
J_A_PMOL_CM2 = 1.0                      # pmol/cm²/impulse (raw reported unit)
J_A = J_A_PMOL_CM2 * 1e-8               # mol/m² (conversion: 1 pmol/cm² = 1e-8 mol/m²)

# ----------------------------------------------------------------------
# A⚠️-type · time-averaged firing rate (microneurography burst structure, no local full text, PMID to be added)
#   ~1 burst/s × 1–3 AP/burst ≈ 2 AP/s (main text [a115][a119])
# ----------------------------------------------------------------------
RATE = 2.0                  # AP/s

# ----------------------------------------------------------------------
# A-type · AP repolarization window (half-width)
#   Belluzzi & Sacchi 1991, Prog Biophys Mol Biol 55:1–30, PMID 2057576:
#   37°C sympathetic postganglionic AP half-width 1.5–2.5 ms → take median 2 ms.
# ----------------------------------------------------------------------
TAU_RELEASE = 2e-3          # s = 2 ms


def varicosity_area(r: float, c: float) -> float:
    """B-type · single varicosity surface area S_v (prolate spheroid formula).

    Raw measured value → formula: the varicosity is nearly a prolate spheroid, equatorial radius r = d_v/2, polar radius c = l_v/2.
      S_v = 2π r² + 2π (r·c/e)·arcsin(e),  e = √(1 − r²/c²)
    Source: geometric model + Goyal & Chaudhury 2013 (PMID 23535140).
    """
    e = math.sqrt(1.0 - (r / c) ** 2)                 # eccentricity (dimensionless)
    term = 2.0 * math.pi * r ** 2 + 2.0 * math.pi * (r * c / e) * math.asin(e)
    return term


def q_per_ap(S_v: float) -> float:
    """B-type · single varicosity per-AP K⁺ efflux Q_b = j_a · S_v.

    Raw measured value → formula: axonal per-AP efflux j_a (Keynes & Ritchie 1965, PMID 5853895)
      × single varicosity surface area S_v (prolate spheroid, Goyal & Chaudhury 2013, PMID 23535140).
    """
    return J_A * S_v                              # [mol/AP]


def instantaneous_peak(j_a: float, a: float, tau_release: float) -> float:
    """B-type · source a instantaneous peak u_inst (within the AP repolarization window, 3D spherical-source surface steady state).

    Raw measured value → formula: u_inst = j_inst · a / D_free, j_inst = j_a/τ_release.
    Parameter sources: j_a (Keynes & Ritchie 1965, PMID 5853895),
      a = 0.5 µm (varicosity radius, Goyal & Chaudhury 2013, PMID 23535140),
      τ_release = 2 ms (AP repolarization duration, standard physiology, to be supplemented locally),
      D_free (Syková & Nicholson 2008, PMID 18923183).
    Note: this peak lasts only ~2 ms and is low-pass filtered out by the terminal membrane τ_m≈110 ms (see design doc §4.4 note),
      and is not used for the "reach threshold" criterion. Returns [mol/m³].
    """
    j_inst = j_a / tau_release                     # [mol/m²/s] instantaneous flux within the repolarization window
    return j_inst * a / D_free


def time_averaged_u(Q_b: float, rate: float, r: float) -> float:
    """B-type · source a time-averaged near-field concentration u_a^ss (3D steady-state point source).

    Raw measured value → formula: u_a^ss = Φ_a / (4π·D_free·r), Φ_a = Q_b·rate.
    Parameter sources: Q_b (§4.2, [2] in this script),
      rate ≈ 2 AP/s (microneurography [a115][a119], to be supplemented locally),
      r = 0.5–2 µm (terminal read distance, sub-document 03),
      D_free (Syková & Nicholson 2008, PMID 18923183).
    Note: clearance length λ_D≈218 µm ≫ r, so the exponential decay term is neglected. Returns [mol/m³], 1 mol/m³ ≡ 1 mM.
    """
    phi_a = Q_b * rate                              # [mol/s] time-averaged source strength
    return phi_a / (4.0 * math.pi * D_free * r)


def compute():
    """Return this module's hand-off quantities (see 01_sympathetic_recruitment.md §5).

    · S_v: single varicosity surface area [m²]
    · Q_b: single varicosity per-AP K⁺ efflux [mol/AP]
    · u_a_ss: source a time-averaged near-field concentration [mol/m³ ≡ mM], keyed by read distance r ∈ {0.5, 1.0, 2.0} µm
    """
    r_med = D_V_MED / 2.0
    c_med = L_V_MED / 2.0
    S_v = varicosity_area(r_med, c_med)
    Q_b = q_per_ap(S_v)
    u_a_ss = {r_um: time_averaged_u(Q_b, RATE, r_um * 1e-6) for r_um in (0.5, 1.0, 2.0)}
    return {"S_v": S_v, "Q_b": Q_b, "u_a_ss": u_a_ss}


if __name__ == "__main__":
    print("=" * 64)
    print("Submodule 01 · Sympathetic recruitment and output — computation output")
    print("=" * 64)

    # ---- B-type · single varicosity surface area (median geometry) ----
    r_med = D_V_MED / 2.0        # 0.5 µm
    c_med = L_V_MED / 2.0        # 1.0 µm
    S_v = varicosity_area(r_med, c_med)
    S_v_min = varicosity_area(D_V_MIN / 2.0, L_V_MIN / 2.0)
    S_v_max = varicosity_area(D_V_MAX / 2.0, L_V_MAX / 2.0)
    print(f"\n[1] single varicosity surface area S_v (prolate spheroid, r={r_med*1e6:.2f} µm, c={c_med*1e6:.2f} µm)")
    print(f"    S_v = {S_v*1e12:.2f} µm²   (range {S_v_min*1e12:.2f}–{S_v_max*1e12:.2f} µm²)")
    print(f"    vs v1 bare constant 1.13 µm²: median ratio = {S_v*1e12/1.13:.2f}×")

    # ---- B-type · single varicosity per-AP K⁺ efflux ----
    Q_b = q_per_ap(S_v)
    print(f"\n[2] single varicosity per-AP K⁺ efflux Q_b = j_a · S_v")
    print(f"    j_a = {J_A_PMOL_CM2} pmol/cm²/impulse (Keynes & Ritchie 1965)")
    print(f"    Q_b = {Q_b:.3e} mol/AP   (vs v1's 1.13e-20, ratio = {Q_b/1.13e-20:.2f}×)")
    print(f"    ≈ Q_b·N_A = {Q_b*N_A:.2e} K⁺ atoms/AP (single-varicosity single-AP atom count, not a free parameter)")

    # ---- B-type · capacitance-independent check ----
    C_M = 1e-2                     # F/m², specific membrane capacitance 1 µF/cm² (standard physiology, no single PMID)
    DV_AP = 0.1                    # V, AP amplitude ~100 mV (standard physiology)
    F = 96485.0                    # C/mol, Faraday constant (SI-defined quantity)
    j_a_check = C_M * DV_AP / F    # mol/m²
    print(f"\n[3] per-AP efflux capacitance-independent check j_a' = C_m·ΔV/F")
    print(f"    j_a' = {j_a_check*1e8:.2f} pmol/cm²  (vs Keynes & Ritchie 1.0, relative difference "
          f"{abs(j_a_check/J_A - 1.0)*100:.1f}%)")

    # ---- B-type · source a instantaneous peak (clarification only, not used for the reach-threshold criterion) ----
    A_V = r_med                    # 0.5 µm, varicosity radius (Goyal & Chaudhury 2013)
    u_inst = instantaneous_peak(J_A, A_V, TAU_RELEASE)
    print(f"\n[4a] source a instantaneous peak u_inst = j_inst·a/D_free (AP repolarization 2 ms window)")
    print(f"     u_inst = {u_inst*1e3:.2f} mM  (only ~2 ms, low-pass filtered out by the terminal membrane τ_m≈110 ms)")

    # ---- B-type · source a time-averaged near-field concentration (for the reach-threshold criterion) ----
    u_a_ss = compute()["u_a_ss"]
    print(f"\n[4b] source a time-averaged near-field concentration u_a^ss = Φ_a/(4π·D_free·r), rate={RATE:.0f} AP/s")
    for r_um in (0.5, 1.0, 2.0):
        u = u_a_ss[r_um]                       # mol/m³ ≡ mM
        print(f"     r = {r_um:4.1f} µm  →  u_a^ss = {u*1e6:5.2f} nM   ({u*1e3:.4f} µM)")
    u_05 = u_a_ss[0.5]
    U_TH = 6.8                                  # mM = 6.8 mol/m³ (module 03 reach threshold; for display comparison, do not change the value)
    print(f"     → time average is nM-level, ~{math.log10(U_TH/u_05):.0f} orders of magnitude from u_th=6.8 mM")

    print("\n" + "=" * 64)
    print("Hand-off: S_v, Q_b, u_a^ss(2–9 nM); N to be anchored (see design doc §6).")
    print("=" * 64)
