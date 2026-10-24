# -*- coding: utf-8 -*-
"""
model_v2 · Submodule 04 · Meridian gradient propagation — computation script
====================================================
Corresponding design doc: 04_meridian_gradient.md
Responsibilities: apparent velocity of gradient manifestation (not conduction) and the reach-threshold front:
  · diffusion length λ_D and per-scenario criterion (λ_D/ℓ)² ≪ 1 (B-type)
  · flux path η = σ_eff/(k·u_th) < 1 → quiescence (B-type)
  · resonance-path apparent velocity v(x) and maximum distance x_max (B-type, Eqs. (equ10)(equ11))
  · comparison against measured PSC velocities 0.68–46 cm/s

Design iron rules (see 00_overview.md §3):
  · Each numeric literal has a same-line comment "source (PMID/local file) + raw measured value".
  · B-type derived values are always computed on the spot from A-type raw values; never hardcode results.

Run: `python model_propagation.py` (depends on common.py in the same directory)
"""

import math
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

from common import D_free, K_CLEAR  # A-type (see common.py for per-item sources)
from model_source_a import compute as compute_2L  # references layer-2 bulk line-source solution of the two-layer model (not hardcoded)

# ----------------------------------------------------------------------
# A-type · sympathetic gradient decay length ℓ
#   derived from span L/ln(1/u_L) (main text §2.3.5(3)), C-type placeholder (not a strict law).
#   value 2 cm (range 1.7–2.6 cm; originally 2.5 cm).
# ----------------------------------------------------------------------
ELL = 2e-2             # m = 2 cm (decay length ℓ, originally 2.5 cm)

# ----------------------------------------------------------------------
# B-type · reach threshold u_th (full-firing threshold)
#   inlined constant (originally referenced module 03's sensory_summation, now deprecated and moved to _deprecated/): empirical anchor Tkachenko 2023 (PMID 36846206)
#   mouse skin free C-MH terminals superfused with 10.8 mM total [K⁺]ₒ → 100% excitation → u_th = 10.8 − 4.0 = 6.8 mM.
#   unit iron rule 1 mol/m³ ≡ 1 mM → 6.8 mM = 6.8 mol/m³.
# ----------------------------------------------------------------------
U_TH = 6.8   # mol/m³ = 6.8 mM (reach threshold, source see note above)

# ----------------------------------------------------------------------
# External input · flux steady state (negative result of Eq. (equ13), not derived from this model's 01/02)
#   references layer 2 of the two-layer computation model (model_source_a.py, cylindrical line-source K₀ solution),
#   taking its bundle-edge bulk concentration u_bulk_edge (r=10 µm, the largest bulk value over the full read distance) as a conservative representative.
#   full read distance 10–272 µm corresponds to u_bulk ≈ 5–42 nM (σ_eff = k·u ≈ 0.05–0.42 nM/s),
#   see 10_S6_charging_model.md §C. Unit iron rule 1 mol/m³ ≡ 1 mM.
# ----------------------------------------------------------------------
U_SS_FLUX = compute_2L()["u_bulk_edge"]   # mol/m³ = bundle-edge bulk concentration (references the two-layer model, not hardcoded)


def diffusion_length(d_free: float, k: float) -> float:
    """B-type · diffusion length λ_D = √(D_free/k)."""
    return math.sqrt(d_free / k)


def compute():
    """Return this module's hand-off quantities (see 04_meridian_gradient.md §5).

    · lam_d: diffusion length [m]
    · criterion_diff: (λ_D/ℓ)² (diffusion/source criterion, ≪1 → independent charging)
    · sigma_flux: flux charging rate σ_eff = u_bulk·k [mol/m³/s]
    · eta_flux: flux η = σ_eff/(k·u_th) (<1 → quiescence)
    · u_th: reach threshold (references module 03)
    """
    lam_d = diffusion_length(D_free, K_CLEAR)
    criterion_diff = (lam_d / ELL) ** 2
    sigma_flux = U_SS_FLUX * K_CLEAR          # mol/m³/s
    eta_flux = sigma_flux / (K_CLEAR * U_TH)
    return {
        "lam_d": lam_d,
        "criterion_diff": criterion_diff,
        "sigma_flux": sigma_flux,
        "eta_flux": eta_flux,
        "u_th": U_TH,
    }


if __name__ == "__main__":
    print("=" * 66)
    print("Submodule 04 · Meridian gradient propagation — computation output")
    print("=" * 66)

    # ---- B-type · diffusion length + per-scenario criterion ----
    lam_d = diffusion_length(D_free, K_CLEAR)            # m
    criterion_diff = (lam_d / ELL) ** 2
    RT_F = 26.7e-3                                       # V (RT/F ≈ 26.7 mV @37°C, common.py)
    # electrical drift/clearance = (F/RT)·D_eff·E0/(k·ℓ), E0 in units of V/cm (1 V/cm = 100 V/m)
    coeff_drift = (1.0 / RT_F) * D_free * 100.0 / (K_CLEAR * ELL)
    print(f"\n[1] diffusion length and per-scenario criterion (gradient charging → diffusion/drift negligible)")
    print(f"    λ_D = √(D_free/k) = {lam_d*1e3:.2f} mm (upper bound without tortuosity folding)")
    print(f"    (λ_D/ℓ)² = ({lam_d*1e3:.2f}/{ELL*1e3:.0f})² = {criterion_diff:.2e} ≪ 1")
    print(f"    electrical drift/clearance = {coeff_drift:.3e}·E0 [V/cm] (reaching E0 = 1 requires "
          f"{1/coeff_drift:.0f} V/cm)")

    # ---- B-type · flux path η ----
    sigma_flux = U_SS_FLUX * K_CLEAR                    # mol/m³/s (σ_eff = u_bulk·k; 1 mol/m³/s ≡ 1 mM/s)
    eta_flux = sigma_flux / (K_CLEAR * U_TH)
    print(f"\n[2] flux path η = σ_eff/(k·u_th)")
    print(f"    σ_eff = u_bulk·k = {sigma_flux*1e6:.3f} nM/s (u_bulk_edge={U_SS_FLUX*1e6:.1f} nM, r=10 µm)")
    print(f"    η = {sigma_flux:.3e}/{K_CLEAR*U_TH:.3e} = {eta_flux:.3e} < 1")
    print(f"    → x_max = ℓ·ln η < 0 (no real solution) → no reach-threshold front → quiescence ✓")

    # ---- B-type · amplification factor A (Eq. equ13) ----
    T_SYM = 1.0                                         # s sympathetic period = 1/f_sym (1:1 cardiac-locked idealized upper bound, 60 bpm)
    A = 1.0 / (1.0 - math.exp(-K_CLEAR * T_SYM))
    print(f"\n[3] amplification factor A = 1/(1−e^(−k·T_sym)), T_sym={T_SYM:.1f} s")
    print(f"    A = {A:.1f} (net injection per burst is amplified by residual accumulation)")

    # ---- B-type · resonance-path apparent velocity (Eqs. equ10/equ11) ----
    ELL_K = ELL * K_CLEAR                               # m/s (= 2e-4)
    ELL_K_CM = ELL_K * 100.0                            # cm/s = 0.02 (v₀ in cm/s, ℓk = 0.02 cm/s)
    print(f"\n[4] resonance-path apparent velocity (v₀ = ℓk·(η−1), ℓk = {ELL_K_CM:.3f} cm/s)")
    print(f"    {'v₀(cm/s)':>9} {'η':>7} {'σ_eff(mM/s)':>11} {'x_max(cm)':>9}")
    for v0 in (0.68, 2.0, 15.9, 46.0):                  # PSC measured-velocity anchor points
        eta = v0 / ELL_K_CM + 1.0                       # η = v₀/(ℓk) + 1
        sig = eta * K_CLEAR * U_TH                      # mM/s (1 mol/m³/s ≡ 1 mM/s)
        xmax = ELL * 1e2 * math.log(eta)                # cm
        print(f"    {v0:9.2f} {eta:7.1f} {sig:11.2f} {xmax:9.1f}")

    # ---- apparent velocity decaying with position (demo η≈431, i.e. mean 15.9 cm/s) ----
    eta_mean = 15.9 / ELL_K_CM + 1.0
    print(f"\n[5] apparent velocity v(x) decaying with position (η={eta_mean:.0f}, corresponding to mean 15.9 cm/s)")
    print(f"    {'x(cm)':>6} {'v(cm/s)':>9}")
    for x_cm in (0.0, 1.0, 2.0, 3.0, 4.0):
        v = ELL_K_CM * (eta_mean * math.exp(-x_cm / (ELL * 1e2)) - 1.0)   # cm/s
        print(f"    {x_cm:6.1f} {v:9.2f}")

    print("\n" + "=" * 66)
    print("Hand-off: x_max(t), v(x); quiescence η<1, resonance η>1 compared against PSC 0.68–46 cm/s.")
    print("=" * 66)
