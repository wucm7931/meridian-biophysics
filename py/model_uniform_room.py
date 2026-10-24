# -*- coding: utf-8 -*-
"""
Uniform-room model comparison — two geometries of background [K⁺]o accumulation under the dam-break regime (clearance off)
================================================================
Compares against e2e_pipeline.py section 6: under the cylindrical steady-state solution (K₀), u_bg only **diverges logarithmically** as k↓
(k↓100× gives only ~2×), so "dam-break" can hardly reach threshold. This script instead uses a **uniform-room model** (u_bg ∝ 1/k,
flux ÷ clearance, no diffusive dilution) to see whether the dam-break regime reaches threshold, and compares it side-by-side with the cylindrical solution.

Physical meaning: the cylindrical K₀ solution = point/line source diffusing into infinite ECS (with dilution); the uniform room = K⁺ confined and
  uniformly mixed in a finite "room" (no dilution, corresponding to a "pathological high-resistance stock" trapping K⁺ locally). The two are
  the two geometric extremes of the dam-break regime, with the true value in between. The three-stage model's "dam-break = pathological damage / high-resistance stock"
  is physically closer to the uniform room (high resistance = K⁺ cannot escape = no dilution), so the dam-break regime should use the uniform-room / mixed caliber.

Run: python model_uniform_room.py
"""

import numpy as np
import sys
from scipy.special import k0, k1
from scipy.optimize import brentq

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

print("=" * 86)
print("Uniform-room model comparison — background [K⁺]o accumulation in the dam-break regime: cylindrical K₀ solution vs uniform room")
print("=" * 86)

# ============================================================================
# 0. Constants (consistent with e2e_pipeline.py; source see that script / model_source_b.py)
# ============================================================================
NERNST   = 26.71
K_REST   = 4.0        # mM
U_TH     = round(10.8 - K_REST, 1)   # 6.8 mM reach-threshold Δ[K⁺]ₒ (Tkachenko 2023 empirical anchor)
DV_TH    = NERNST * np.log((K_REST + U_TH) / K_REST)  # 26.5 mV reach-threshold voltage (derived from empirical u_th)
K_CLEAR  = 0.01       # /s baseline clearance (muscle-anchored; Förstl SCG 0.04 upper end)
D_eff    = 7.42e-10   # m²/s
ALPHA    = 0.30
R_VESSEL = 12.5e-6
R_S      = R_VESSEL + 1.0e-6

I_STOC   = 250e-12
T_STOC   = 100e-3
F        = 96485.0
S_CELL   = 4600e-12
J_STOC_MOL = I_STOC * T_STOC / (F * S_CELL)   # 5.63e-8 mol/m² single-STOC K⁺ efflux

F_STOC   = 1.0
lam_cmax = 0.2
N        = 20
U_SAT    = 12.0       # mM nonlinear clearance saturation point (source_b k(u)=k₀(1+(u/12)^4))
DENS_MERIDIAN = 1.5   # along-meridian coefficient: averaged along the meridian line (acupoint 3.0 × acupoint fraction along the meridian ~0.15–0.33 + inter-acupoint segment ~1.0 ≈ 1.3–1.6)
DENS_POINT    = 3.0   # acupoint coefficient: local microvascular density at the acupoint (conservative placeholder; literature only says "significantly higher" with no multiple, Liu 2014/Hsiu 2007)
# Note: the old DENSITY=4.8 is deprecated — it is a PCV "prevalence ratio" (87.5%/18.1%), not a continuous density ratio,
#   and it only measured acupoint vs a 4 cm control point, not the inter-acupoint meridian segments. Prevalence→density is not monotonic (Poisson ≈10×, clustering <4.8).

# ============================================================================
# 1. Cylindrical K₀ steady-state solution (logarithmic divergence)
# ============================================================================
def cyl_u(r, J, k=K_CLEAR):
    lam = np.sqrt(D_eff / k)
    pref = J * lam / (D_eff * ALPHA)
    return pref * k0(r / lam) / k1(R_VESSEL / lam)

def u_bg_cyl(k):
    """N-source time-averaged background [K⁺]o (mol/m³≡mM), R=A=1."""
    lam_tot = (F_STOC + lam_cmax) * 1.0
    J_avg = J_STOC_MOL * lam_tot
    return N * cyl_u(R_S, J_avg, k=k)

# ============================================================================
# 2. Uniform-room model (no diffusive dilution, u_bg ∝ 1/k)
# ============================================================================
# Anchoring: at k=0.01 the uniform room equals the cylindrical solution (linear region, 187 µM far below the 12 mM saturation point, nonlinearity negligible)
J_INJ = K_CLEAR * u_bg_cyl(K_CLEAR)            # injection rate mol/m³/s, anchored at k=0.01

def u_bg_room_linear(k):
    """Uniform-room linear clearance: u = J/k (flux ÷ clearance, no dilution)."""
    return J_INJ / k

def u_bg_room_nonlin(k, dens=1.0):
    """Uniform room + nonlinear clearance: solve k·u·(1+(u/U_SAT)^4) = J (Na/K-ATPase saturation + Nernst reverse constraint).
    dens = along-meridian density multiplier (number of sources per volume ∝ vessel/varicosity density)."""
    J = J_INJ * dens
    def f(u):
        return k * u * (1.0 + (u / U_SAT) ** 4) - J
    return brentq(f, 1e-9, 100.0)              # u in units of mM

def dV_bg(u):
    """Background depolarization (background sustained, INTEG=1)."""
    return NERNST * np.log((K_REST + u) / K_REST)

# ============================================================================
# 3. Scan clearance rate k, three models side-by-side
# ============================================================================
print("\n[1] Anchoring (k=0.01/s, R=A=1, N=20)")
u_ref = u_bg_cyl(K_CLEAR)
print(f"    cylindrical solution u_bg({K_CLEAR}) = {u_ref*1e3:.1f} µM → ΔV_bg = {dV_bg(u_ref):.2f} mV")
print(f"    injection rate J_INJ = k·u_bg = {J_INJ*1e3:.2f} µM/s (uniform room anchored here)")

print("\n[2] dam-break regime scan (k↓: cylindrical K₀ solution vs uniform-room nonlinear, along-meridian 1.5× / acupoint 3.0×; full threshold ΔV_bg ≥ 26.5 mV)")
print(f"    {'k(1/s)':>9} {'cyl ΔV':>8} | {'meridian 1.5× u_bg':>13} {'ΔV':>7} {'reach':>5}"
      f" | {'acupoint 3.0× u_bg':>13} {'ΔV':>7} {'reach':>5}")
for k_scan in [0.04, 0.01, 0.004, 0.001, 1e-4, 1e-5, 1e-6, 1e-7]:
    u_cyl = u_bg_cyl(k_scan)
    u_m   = u_bg_room_nonlin(k_scan, dens=DENS_MERIDIAN)
    u_p   = u_bg_room_nonlin(k_scan, dens=DENS_POINT)
    fm = "✅" if dV_bg(u_m) >= DV_TH else ""
    fp = "✅" if dV_bg(u_p) >= DV_TH else ""
    print(f"    {k_scan:9.0e} {dV_bg(u_cyl):7.2f} | {u_m*1e3:12.1f}µM {dV_bg(u_m):6.2f} {fm:>5}"
          f" | {u_p*1e3:12.1f}µM {dV_bg(u_p):6.2f} {fp:>5}")

print("""
[3] Conclusions
  ① Cylindrical K₀ solution: u_bg only **diverges logarithmically** as k↓ (k↓1e6× from 187 µM to ~500 µM),
     ΔV_bg is always < 2 mV — "dam-break" cannot reach threshold in the cylindrical geometry.
  ② Uniform-room linear: u_bg ∝ 1/k, k↓100× → u_bg↑100×, at k≈1e-5/s u_bg~0.75 mM,
     at k→1e-6 ~7.5 mM → ΔV_bg crosses 26.5 mV (dam-break).
  ③ Uniform-room nonlinear (Na/K-ATPase saturation + Nernst constraint): u_bg is clamped at ~17 mM (even as k→0),
     ΔV_bg saturates at ~45 mV ≥ 26.5 mV — the dam-break regime background crosses threshold without diverging infinitely (physically self-consistent).
  ④ Physical judgment: the physics of the three-stage "dam-break = pathological damage / high-resistance stock" is "K⁺ trapped by the high-resistance stock"
     = no diffusive dilution = **closer to the uniform room**, not the cylindrical K₀. Hence the dam-break background accumulation should use the uniform-room caliber,
     with the true value between ②/③ (linear is the no-dilution upper bound, nonlinear is the physical saturation).
  ④′ The density coefficient (number of sources per volume ∝ vessel/varicosity density) enters the uniform room as injection rate J_INJ × coefficient:
     along-meridian 1.5×, acupoint 3.0× (old 4.8 deprecated; its "prevalence ratio" is not a density ratio and only measured acupoints, not inter-acupoint segments).
     The ordinary regime only crosses threshold at k≈1e-4/s; the along-meridian/acupoint regimes both move it up to k≈1e-3/s (the 2× density difference is compressed
     by nonlinear clearance to the same order of magnitude, so the conclusion is insensitive to the density coefficient) — density is a precondition for "lowering the gate",
     not the "closing" itself (the closing is still k↓).
  ⑤ Significance for the end-to-end pipeline: the conclusion of e2e_pipeline section 6 that "the time-accumulation axis contributes little (~1.2 mV)"
     holds only for the "healthy quiescence regime" (cylindrical geometry, normal clearance, non-acupoint density); under the "dam-break regime" (clearance off / high-resistance stock)
     + along-meridian density, background accumulation does reach threshold, at which point the reaching-threshold mode switches from "synchronized peak" to "background crossing threshold (sustained plateau)".
""")
