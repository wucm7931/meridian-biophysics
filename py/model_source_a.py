# -*- coding: utf-8 -*-
"""
Source a · bulk line-source steady state — §S6 §C "two-layer computation model" layer 2 (negative result)
=====================================================================

【Position of this script in the model】
  Corresponds to 10_S6_charging_model.md §C "flux cannot reach threshold (negative result)" layer 2:
  bottom-up reconstruction of the steady-state concentration u_bulk in the bulk extracellular space (ECS) produced by K⁺ efflux from sympathetic varicosities in the "quiescent state" (healthy resting, non-pathological/non-damaged/non-acupoint-along-meridian), proving it cannot reach
  the PSC threshold u_th = 6.8 mM — i.e. the quantitative answer for "no PSC in a healthy person" (negative result = quiescence baseline).

【Two-layer computation model (aligned with §S6 §C)】
  Layer 1 (soma/channel level):
    BSW 1991 five-conductance Hodgkin-Huxley model, checking the channel model's RMP / cleft
    Δ[K⁺] / per-AP charge, only for channel validation, not as a single-varicosity source term.
  Layer 2 (this script, population bulk line source) — six bottom-up steps:
    [1] Q_b  = j_a × S_v                    single-varicosity per-AP K⁺ (axon efflux anchored)
    [2] n_b  = varicosity density × fiber count             varicosity line density (Luff 1987 × fibers per bundle)
    [3] rate = 0.33 Hz (measured average firing rate)   resting single-fiber firing rate (Macefield & Wallin 1999)
    [4] Φ    = Q_b × n_b × rate            line-source strength (K⁺ mol / per meter / per second)
    [5] u(r) = Φ/(2π D_eff) × K₀(r/λ_D)    cylindrical line-source reaction-diffusion steady state (K₀ solution)
    [6] σ_eff= k × u(r)                    equivalent net charging rate (steady state u(r) = σ_eff/k)

【Why a cylindrical line source (K₀ solution) and not a 0D well-mixed】
  A sympathetic bundle is a 1-D line running along an arteriole (~40 fibers, radius ~15 µm). The K⁺ released by a single varicosity
  diffuses two-dimensionally across the cross-section while being cleared at first order (rate k). The steady state of this "line source + diffusion
  + clearance" is not the 0D approximation that uniformly spreads the source over some radius (the latter differs by 1–2 orders of magnitude from the arbitrary dilution
  radius), but has the closed-form solution u(r) = Φ/(2πD_eff)·K₀(r/λ_D)——
  where λ_D = √(D_eff/k) is the diffusion-clearance length and K₀ is the zero-order modified Bessel function of the second kind.

【Design iron rules (00_overview.md §3)】
  · Every numeric literal has an inline comment "provenance (PMID/local file) + raw measurement value".
  · All B-class derived values are computed on the spot from A-class raw values; results must not be hardcoded.
  · Unit iron rule: 1 mol/m³ ≡ 1 mM; 1 pmol/cm² = 1e-8 mol/m²; 1 µm² = 1e-12 m².

Run: `python model_source_a.py` (depends on common.py in the same directory)
"""

import math
import sys

# Windows console defaults to GBK and cannot encode µ/²/× etc.; force stdout to UTF-8 output.
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

# A/B class · shared tissue parameters (single source of truth, see common.py comments)
#   D_eff = D_free/λ² = 7.42×10⁻¹⁰ m²/s (D_free divided by tortuosity, physically correct value);
#   K_CLEAR = k = 0.01 s⁻¹ (C-class main value, range [0.003, 0.04], muscle/spinal-cord anchored).
from common import D_eff, K_CLEAR, N_A, NERNST, K_REST

# ----------------------------------------------------------------------
# A class · axon per-AP K⁺ efflux (raw measurement value + trunk→varicosity extrapolation assumption)
#   Keynes & Ritchie 1965, J Physiol 179:333–367, PMID 5853895:
#   mammalian unmyelinated nerve trunk (rabbit vagus C fiber) per-impulse K⁺ efflux = 1.0 pmol/cm².
#   ⚠️ This value was measured on a nerve trunk; carrying it over to the varicosity membrane (same axon, same membrane, same AP, same per-area efflux)
#      is a first-order assumption, not independently verified on varicosity swellings; the upper bound is locked by charge conservation C_m·ΔV/F≈1.04 pmol/cm².
# ----------------------------------------------------------------------
J_A_PMOL_CM2 = 1.0               # pmol/cm²/impulse (raw reported unit)
J_A = J_A_PMOL_CM2 * 1e-8        # mol/m² (conversion: 1 pmol/cm² = 1e-8 mol/m²)

# ----------------------------------------------------------------------
# A class · varicosity geometry (raw measurement values)
#   Goyal & Chaudhury 2013, Auton Neurosci 176:11–31, PMID 23535140:
#   sympathetic varicosities are spindle-shaped, diameter d_v = 0.5–2.0 µm, length l_v = 1–3 µm (take median).
# ----------------------------------------------------------------------
D_V_MED = 1.0e-6                 # median diameter 1.0 µm [m]
L_V_MED = 2.0e-6                 # median length 2.0 µm [m]

# ----------------------------------------------------------------------
# A class · resting sympathetic single-fiber average firing rate (raw measurement value, cited directly)
#   Macefield & Wallin 1999, J Physiol 516:293–301, PMID 10066942 (main text [c97]):
#   resting (high MSNA group) single-unit average firing rate = 0.33 ± 0.04 Hz, i.e. about 0.33 AP per fiber per second;
#   structure is 35% cardiac-cycle participation, 78% of firing intervals are single spikes (s_avg ≈ 1).
#   ⚠️ This rate is a "fiber firing rate", not a vesicle release probability — the K⁺ repolarization efflux happens on every
#      AP that invades the varicosity, with no vesicle-probability gating (design iron rule, see sub-design 01 §2).
#   (The old Macefield 1994, PMID 7707244 p_fire=0.21 has been replaced by the 1999 version,
#     see 03_hypotheses_grid_psc_ref.md [c97]-backup.)
# ----------------------------------------------------------------------
RATE = 0.33                      # Hz = AP/s (resting single-fiber average firing rate, directly measured)

# ----------------------------------------------------------------------
# C class · varicosity line density (placeholder/to be anchored)
#   Luff SE, McLachlan EM, Hirst GDS 1987, J Comp Neurol 257:578–594,
#     PMID 3693598: sympathetic varicosity density 7.6–9.5 per 100µm/fiber (take 8.0).
#   The fiber count per arteriolar sympathetic bundle N_FIBERS ≈ 40 is a model placeholder (C class, to be anchored:
#   "how many varicosities each sensory terminal faces" needs serial ultrathin-section microscopic counting; no precise value locally).
#   → n_b = 8 varicosities/100µm/fiber × 40 fibers = 3.2×10⁶ varicosities/m.
# ----------------------------------------------------------------------
N_VAR_100UM = 8.0                # varicosities/100µm/fiber (Luff 1987, median of 7.6–9.5)
N_FIBERS = 40.0                  # fibers per bundle (C-class placeholder, to be anchored)

# ----------------------------------------------------------------------
# B class · sensory-terminal activation threshold u_th (reach-threshold Δ[K⁺]ₒ, computed on the spot, not hardcoded)
#   Empirical anchor (Tkachenko 2023, PMID 36846206): superfusion of the receptive field of mouse skin free C-MH terminals with
#     10.8 mM total [K⁺]ₒ → 100% excitation (23/23) → u_th = total threshold − K_REST = 10.8 − 4.0 = 6.8 mM.
#   (The old Nernst back-solution 30.4 mV → 8.5 mM is dropped: DRG distal axon + K⁺ selectivity assumption, overestimating ~25%.)
#   Unit iron rule: 1 mol/m³ ≡ 1 mM → 6.8 mM = 6.8 mol/m³.
# ----------------------------------------------------------------------
U_TH = round(10.8 - K_REST, 1)                    # mol/m³ = 6.8 mM


def varicosity_area(r: float, c: float) -> float:
    """B class · single-varicosity surface area S_v (prolate spheroid formula).

    Raw measurement value → formula: the varicosity approximates a prolate spheroid, equatorial radius r = d_v/2, polar radius c = l_v/2.
      S_v = 2π r² + 2π (r·c/e)·arcsin(e),  e = √(1 − r²/c²)
    Provenance: geometric model + Goyal & Chaudhury 2013 (PMID 23535140).
    """
    e = math.sqrt(1.0 - (r / c) ** 2)                 # eccentricity (dimensionless)
    return 2.0 * math.pi * r ** 2 + 2.0 * math.pi * (r * c / e) * math.asin(e)


def k0(x: float) -> float:
    """K₀(x): zero-order modified Bessel function of the second kind (the core special function of the cylindrical line-source steady-state solution).

    Appears in the "2D line source + first-order clearance" steady state u(r) = Φ/(2πD)·K₀(r/λ).
    Uses the Abramowitz & Stegun 1964 §9.8 polynomial approximation (|ε| < 1.5×10⁻⁷):
      0 < x ≤ 2:  K₀(x) = −ln(x/2)·I₀(x) − γ + Σ a_k (x/2)^{2k}
      x ≥ 2:      √x·eˣ·K₀(x) = Σ b_k (2/x)^k
    Parameters are only read in the range r/λ_D ≲ 2, where the main branch (x≤2) suffices; the x>2 branch is for extrapolation.
    """
    if x <= 2.0:
        # I₀(x) (zero-order modified Bessel of the first kind), A&S 9.8.1
        t = x / 3.75
        t2 = t * t
        i0 = (1.0 + 3.5156229 * t2 + 3.0899424 * t2 ** 2 + 1.2067492 * t2 ** 3
              + 0.2659732 * t2 ** 4 + 0.0360768 * t2 ** 5 + 0.0045813 * t2 ** 6)
        y = x / 2.0
        y2 = y * y
        # K₀(x), A&S 9.8.5 (γ = 0.57721566)
        return (-math.log(y) * i0
                + (-0.57721566 + 0.42278420 * y2 + 0.23069756 * y2 ** 2
                   + 0.03488590 * y2 ** 3 + 0.00262698 * y2 ** 4
                   + 0.00010750 * y2 ** 5 + 0.00000740 * y2 ** 6))
    else:
        # K₀(x), A&S 9.8.6
        t = 2.0 / x
        poly = (1.25331414 - 0.07832358 * t + 0.02189568 * t ** 2
                - 0.01062446 * t ** 3 + 0.00587872 * t ** 4
                - 0.00251540 * t ** 5 + 0.00053208 * t ** 6)
        return math.exp(-x) / math.sqrt(x) * poly


def bulk_line_source_profile():
    """Return the per-item intermediate quantities of the cylindrical line-source steady-state solution u(r) = Φ/(2πD_eff)·K₀(r/λ_D).

    Returns a dict:
      S_v       single-varicosity surface area [m²]
      Q_b       single-varicosity per-AP K⁺ [mol/AP]
      n_b       varicosity line density [/m]
      rate      resting firing rate [AP/s]
      Phi       line-source strength [mol/(m·s)]
      lambda_D  diffusion-clearance length [m]
      prefactor Φ/(2πD_eff) [mol/m³ ≡ mM]
    """
    r_med = D_V_MED / 2.0
    c_med = L_V_MED / 2.0
    S_v = varicosity_area(r_med, c_med)              # [1] single-varicosity area
    Q_b = J_A * S_v                                  # single-varicosity per-AP K⁺
    n_b = N_VAR_100UM / 100e-6 * N_FIBERS            # [2] varicosity line density
    rate = RATE                                      # [3] resting firing rate (directly measured)
    Phi = Q_b * n_b * rate                           # [4] line-source strength
    lambda_D = math.sqrt(D_eff / K_CLEAR)            # diffusion-clearance length
    prefactor = Phi / (2.0 * math.pi * D_eff)        # concentration scale (coefficient in front of K₀)
    return {"S_v": S_v, "Q_b": Q_b, "n_b": n_b, "rate": rate,
            "Phi": Phi, "lambda_D": lambda_D, "prefactor": prefactor}


def compute():
    """Module hand-off quantities (referenced by 04_meridian_gradient, not hardcoded).

    Return the bulk concentration and equivalent net charging rate of the cylindrical line-source steady-state solution u(r) = Φ/(2πD_eff)·K₀(r/λ_D) at the out-of-bundle read distances
    r ∈ {10, 30, 90, 272} µm, plus the summary value at the "bundle edge" (r=10 µm, the largest bulk value closest to the bundle).

    · u_bulk: {read distance µm: concentration [mol/m³ ≡ mM]}
    · sigma_eff: {read distance µm: k·u [mol/m³/s]}
    · u_bulk_edge: concentration at bundle edge r=10 µm [mol/m³] (conservative representative value — if even the largest
      bulk value is ≪ threshold, the negative result is robust across all read distances)
    · prefactor / lambda_D: the two scale parameters of the full solution
    """
    p = bulk_line_source_profile()
    read_um = (10.0, 30.0, 90.0, 272.0)
    u_bulk = {}
    sigma_eff = {}
    for r_um in read_um:
        u = p["prefactor"] * k0(r_um * 1e-6 / p["lambda_D"])
        u_bulk[r_um] = u
        sigma_eff[r_um] = K_CLEAR * u
    return {
        "u_bulk": u_bulk,
        "sigma_eff": sigma_eff,
        "u_bulk_edge": u_bulk[10.0],     # bundle edge r=10 µm (largest bulk value)
        "prefactor": p["prefactor"],
        "lambda_D": p["lambda_D"],
    }


if __name__ == "__main__":
    print("=" * 72)
    print("Source a · bulk line-source steady state (§S6 §C two-layer computation model layer 2 · negative result)")
    print("=" * 72)

    p = bulk_line_source_profile()

    # ---- [1] single-varicosity area and per-AP K⁺ (layer-2 anchor, consistent with sub-design 01 §4) ----
    print("\n[1] single-varicosity per-AP K⁺ efflux Q_b = j_a × S_v")
    print(f"    S_v = {p['S_v'] * 1e12:.2f} µm²  (prolate spheroid, d=1.0 / l=2.0 µm, Goyal 2013)")
    print(f"    j_a = {J_A_PMOL_CM2} pmol/cm²/impulse  (Keynes & Ritchie 1965)")
    print(f"    Q_b = {p['Q_b']:.3e} mol/AP  (≈ {p['Q_b'] * N_A:.2e} K⁺/AP)")

    # ---- [2] varicosity line density ----
    print("\n[2] varicosity line density n_b = 8/100µm/fiber × 40 fibers")
    print(f"    n_b = {p['n_b']:.2e} varicosities/m  (Luff 1987 × N_FIBERS=40 placeholder)")

    # ---- [3] resting firing rate ----
    print("\n[3] resting firing rate rate (directly measured, Macefield & Wallin 1999 [c97])")
    print(f"    rate = {RATE:.2f} Hz  (0.33 ± 0.04 Hz; 35% cardiac participation, 78% single spike)")

    # ---- [4] line-source strength ----
    print("\n[4] line-source strength Φ = Q_b × n_b × rate")
    print(f"    Φ = {p['Phi']:.3e} mol/(m·s)")

    # ---- [5] cylindrical line-source steady state (K₀ solution) ----
    print("\n[5] cylindrical line-source reaction-diffusion steady state u(r) = Φ/(2πD_eff) × K₀(r/λ_D)")
    print(f"    D_eff = {D_eff:.3e} m²/s  (D_free/λ², divided by tortuosity)")
    print(f"    λ_D   = √(D_eff/k) = {p['lambda_D'] * 1e6:.0f} µm  (k = {K_CLEAR} s⁻¹)")
    print(f"    Φ/(2πD_eff) = {p['prefactor'] * 1e6:.2f} nM  (concentration scale)")
    print("\n    read distance r        K₀(r/λ_D)     u_bulk(r)      σ_eff = k·u_bulk")
    read_points = (10.0, 30.0, 90.0, 272.0)   # µm: out-of-bundle proximal → clearance length λ_D
    for r_um in read_points:
        r = r_um * 1e-6
        kk = k0(r / p["lambda_D"])
        u = p["prefactor"] * kk              # mol/m³ ≡ mM
        sigma = K_CLEAR * u                  # mol/m³/s = equivalent net charging rate
        print(f"    r = {r_um:5.0f} µm   {kk:10.3f}    {u * 1e6:8.2f} nM      "
              f"{sigma * 1e6:6.3f} nM/s")

    # ---- [6] negative-result determination ----
    print("\n[6] negative-result determination: u_bulk ≪ u_th = 6.8 mM")
    u_min = p["prefactor"] * k0(272e-6 / p["lambda_D"])
    u_max = p["prefactor"] * k0(10e-6 / p["lambda_D"])
    orders_min = math.log10(U_TH / u_max)
    orders_max = math.log10(U_TH / u_min)
    print(f"    u_bulk ≈ {u_min * 1e6:.1f}–{u_max * 1e6:.1f} nM  "
          f"(read distance {read_points[0]:.0f}–{read_points[-1]:.0f} µm)")
    print(f"    ~{orders_min:.1f}–{orders_max:.1f} orders of magnitude below threshold → flux cannot reach threshold")
    print(f"    (k·u_th = {K_CLEAR * U_TH * 1e3:.2f} µM/s, "
          f"σ_eff ≈ {K_CLEAR * u_min * 1e6:.3f}–{K_CLEAR * u_max * 1e6:.3f} nM/s, "
          f"differing by about {math.log10((K_CLEAR * U_TH) / (K_CLEAR * u_max)):.1f}"
          f"–{math.log10((K_CLEAR * U_TH) / (K_CLEAR * u_min)):.1f} orders of magnitude)")

    # ---- robustness note ----
    print("\n" + "-" * 72)
    print("Robustness note (conclusion independent of free parameters)")
    print("  · u_bulk ∝ n_b ∝ N_FIBERS: even if the fiber count N_FIBERS is taken as 40×10=400 (placeholder scaled")
    print("    up another 10×), u_bulk ≈ 51–416 nM, still about 4.2–5.1 orders of magnitude below 6.8 mM.")
    print("  · k range [0.003, 0.04]: λ_D increases as k decreases (slower clearance → farther diffusion → more dilute),")
    print("    so the negative result is robust across the whole range (see memo_k_clearance_muscle.md).")
    print("  · D_eff uses the divided-by-tortuosity value (physically correct); the undivided D_free upper bound gives an even smaller concentration, making the conclusion more robust.")
    print("=" * 72)
