# -*- coding: utf-8 -*-
"""
Near-field K⁺ concentration: instantaneous peak vs steady state — re-derived from the literature
============================================================

Goal: determine whether the K⁺ efflux of a single sympathetic varicosity action potential (AP)
      can depolarize the adjacent sensory C-fiber terminal (gap 20–150 nm) to its
      threshold voltage (~26.5 mV, derived from the empirical u_th).

The core question addressed here (a correction to S6):
    "Is a steady state required to activate the sensory nerve?"
    — Answer: no. The sensory terminal reads the **single-AP instantaneous K⁺ peak**, not
    the 1 Hz continuous-discharge time average (bulk steady state). The two differ by
    ~3 orders of magnitude, but both are far below the firing threshold.

    S6's original "terminal near field ~1–7 nM" is actually the **bulk steady state of a
    cylindrical line source** (time average, 1 Hz continuous discharge), not the near-field
    concentration the terminal reads at the instant of a single AP. The near-field
    instantaneous peak should be ~1–9 µM (see (A)). The conclusion (the flux cannot reach
    threshold) is unchanged, but the magnitude and wording must be corrected from "nM" to
    "µM (instantaneous)".

Three concentrations are computed (each converted to depolarization via Nernst, then compared
with the 26.5 mV threshold):
  (A) Single-AP instantaneous peak — crack quasi-steady state (with lateral escape) — ~µM  ← the core correction here
  (B) Single-AP instantaneous peak upper bound — instantaneous constraint (no escape) — ~0.1–0.5 mM (non-physical upper bound)
  (C) 1 Hz continuous-discharge bulk steady state — cylindrical line source — ~nM  ← S6's old value (mislabeled as near field)

All formulas, parameters, and units are traced to the literature one by one (PMIDs in the
notes of each section).
Run: python calc_nearfield.py
"""

import numpy as np
import sys
from scipy.special import k0

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

print("=" * 84)
print("Near-field K⁺: instantaneous peak vs steady state (re-derived from the literature)")
print("=" * 84)

# ============================================================================
# 0. Physical constants (standard biophysics textbook, Hille, Ion Channels of Excitable Membranes)
# ============================================================================
F      = 96485.0       # Faraday constant [C/mol]
R_gas  = 8.314         # gas constant [J/(mol·K)]
T      = 310.0         # body temperature 37°C [K]
# Nernst slope: ΔE_K = RT/F · ln([K+]new/[K+]rest)
# R·T/F = 8.314*310/96485 = 0.02671 V = 26.71 mV
NERNST = R_gas * T / F * 1000.0      # [mV] = 26.71

C_m    = 1e-6         # specific membrane capacitance [F/cm²] = 1 µF/cm² (approximately constant for all excitable membranes)
dV_AP  = 0.100        # action potential amplitude [V] = 100 mV

# --- Unit conversion (cgs ↔ molar concentration) ---
# 1 mol/cm³ = 1e3 mol/L = 1e3 M = 1e6 mM = 1e9 µM = 1e12 nM
MOL_CM3_TO_mM = 1e6
MOL_CM3_TO_uM = 1e9
MOL_CM3_TO_nM = 1e12

# ============================================================================
# 1. K⁺ efflux j [mol/cm²/impulse]
# ============================================================================
# Formula (charge conservation, universal): per-AP Na⁺ influx = K⁺ efflux (net charge balance)
#     j = C_m · ΔV_AP / F
#   Value: 1e-6 F/cm² × 0.1 V / 96485 C/mol = 1.036e-12 mol/cm² = 1.04 pmol/cm²
# Independent literature measurement: Keynes & Ritchie 1965, J Physiol 179:333–367, PMID 5853895
#     rabbit C fiber 42K tracer K⁺ efflux ≈ 1.0 pmol/cm²/impulse.
# The two agree: charge conservation is the universal upper bound, and the measurement is an independent confirmation.
j = C_m * dV_AP / F            # [mol/cm² per impulse]

# ============================================================================
# 2. Varicosity geometry [cm] — source Goyal & Chaudhury 2013, PMID 23535140 [c142]
# ============================================================================
# Varicosity: diameter 0.5–2.0 µm, length 1–3 µm, spacing 2–10 µm, >25,000 varicosities per axon.
# Taking mid-values: diameter 1.0 µm (semi-minor axis b_v = 0.5 µm), length 2.0 µm (semi-major axis a_v = 1.0 µm),
# treated as a prolate spheroid.
a_v = 1.0e-4                 # semi-major axis [cm] = 1.0 µm
b_v = 0.5e-4                 # semi-minor axis [cm] = 0.5 µm
# Prolate spheroid surface area: S = 2πb² + 2π(a b/e)·arcsin(e), e = sqrt(1 - b²/a²)
e   = np.sqrt(1.0 - (b_v / a_v) ** 2)
S_v = 2 * np.pi * b_v**2 + 2 * np.pi * (a_v * b_v / e) * np.arcsin(e)   # [cm²]
# Value ≈ 5.4e-8 cm² = 5.4 µm²
Q_b = j * S_v                 # total K⁺ efflux per varicosity per AP [mol/impulse] ≈ 5.6e-20 mol

# ============================================================================
# 3. Gap geometry g and apposition-patch radius R_app [cm]
# ============================================================================
# Sympathetic varicosity ↔ sensory terminal gap (three independent sources):
#   50–1000 nm  Aalkjær et al. 2021, Physiol Rev 101:495, PMID 33270533 [c49]
#   20–150 nm   Bennett 1996, Prog Neurobiol 50:505, PMID 9015825        [c143]
#   50–100 nm   Luff 1996, Anat Embryol 193:515, PMID 8737808            [c139] (effector side)
# Scan over 20/50/100/150 nm.
gaps_nm = np.array([20.0, 50.0, 100.0, 150.0])   # [nm]
gaps_cm = gaps_nm * 1e-7                          # [cm]
R_app   = 0.5e-4                                  # apposition-patch radius [cm] = 0.5 µm (≈ b_v)

# ============================================================================
# 4. Diffusion coefficient D and clearance rate k
# ============================================================================
# D_free = 1.9e-5 cm²/s: K⁺ free-water diffusion coefficient (Syková & Nicholson 2008,
#     Physiol Rev, PMID 18923183, the cited D_KCl).
# ECS tortuosity λ = 1.6, volume fraction α = 0.2 (Nicholson & Syková 1998,
#     Trends Neurosci 21:207, PMID 9610885 [c99]); during ischemia α=0.04, λ=2.2.
#   Effective diffusion D_tort = D_free / λ² = 1.9e-5/2.56 = 7.4e-6 cm²/s
D_free = 1.9e-5              # [cm²/s]
lam_T  = 1.6                 # tortuosity λ (dimensionless)
D_tort = D_free / lam_T**2   # [cm²/s] effective diffusion after tortuosity correction
# First-order clearance −k u (lumped: glial uptake + Na/K pump + blood flow):
#   Hoppe 1991, Brain Res 552:106, PMID 1913170 [c112]: t1/2 = 50–100 s
#       → k = ln2/t1/2 = 0.007–0.014 /s (direct microelectrode measurement, diffusion-dominated)
#   Förstl 1982, Neuroscience 7:3221, PMID 7162635 [c88]: SCG (S6 takes τ≈25 s → 0.04/s)
#   Note: k only affects the bulk steady state (C), not the instantaneous peak (A/B) (instantaneous ~ms ≪ clearance ~25–100 s).
k_fast = 0.04                        # Förstl 1982 (SCG, value used by S6) [1/s]
k_slow = np.log(2.0) / 100.0         # Hoppe 1991 lower bound (t1/2=100 s) [1/s]

# ============================================================================
# 5. Sensory terminal electrogenesis parameters
# ============================================================================
# Type-A reference (transduction-layer electrogenesis, not the threshold basis): Vasylyev & Waxman 2012,
#     J Neurophysiol 108:729, PMID 22572942 [c134] — free C-fiber terminal RMP −60.2 mV,
#     AP threshold −29.8 mV, ΔV=30.4 mV (Nav1.7 amplifies subthreshold / Nav1.8 fires the AP).
V_rest = -60.2           # [mV] terminal RMP (Type-A reference)
V_thr  = -29.8           # [mV] AP threshold (Type-A reference)
# Nav1.7 subthreshold-amplification onset ~5 mV ([c134] abstract "from ~−55 mV"; [c137] Mulcahy 2019 threshold
#   channel, ramp current amplifies small depolarizations). This is a lower threshold for the "minimal detectable depolarization":
#   5 mV still requires Δ[K⁺]o ≈ 0.8 mM (Nernst), still far above the ~µM instantaneous peak.
dV_Nav17 = 5.0            # [mV] minimum depolarization at which Nav1.7 begins amplifying subthreshold depolarization
# Membrane time constant τ_m = 110 ms (Bostock 2003, J Physiol 553:649, PMID 12963801 [c135])
tau_m    = 110e-3        # [s]
# K⁺ transient duration τ_pulse ≈ repolarization time (Kv channel open ~2 ms)
tau_pulse = 2e-3         # [s]
# Resting [K⁺]o = 4 mM (mid-value of interstitial fluid 3.5–5 mM [c321]; note Förstl 1982 measured 6 mM in SCG neuropil)
# Key background: baseline is already ~4 mM, so any physiologically meaningful [K⁺]o fluctuation must start at the µM level,
#   "nM" can only be a time average / single-fiber extreme underestimate, and cannot be the near-field instantaneous value (user's intuition, confirmed).
K_rest_mM = 4.0
# Threshold voltage ΔV_th (Type B, derived from the empirical u_th, not V_thr−V_rest):
#   Empirical anchor Tkachenko 2023 (Front Cell Neurosci 17:1131643, PMID 36846206):
#   superfusing 10.8 mM total [K⁺]ₒ over the receptive field of mouse skin free C-MH terminals → 100% excitation →
#   u_th = 10.8 − 4.0 = 6.8 mM; ΔV_th = NERNST·ln((K_rest+u_th)/K_rest) ≈ 26.5 mV.
u_th_mM = 6.8            # [mM] threshold Δ[K⁺]ₒ (empirical anchor)
dV_th   = NERNST * np.log((K_rest_mM + u_th_mM) / K_rest_mM)   # ≈ 26.5 mV

# ============================================================================
# Helper functions
# ============================================================================
def nernst_mV(dK_mM):
    """K⁺ increment dK_mM [mM] → depolarization ΔE_K [mV] (Nernst, 37°C)"""
    return NERNST * np.log((K_rest_mM + dK_mM) / K_rest_mM)

def integrator_filter(dV_mV, tau_pulse_s, tau_m_s):
    """Membrane integrator filter: rectangular pulse of amplitude dV and duration τ_pulse → peak displacement ΔV_peak
    Formula (RC membrane, Hille 2001): ΔV_peak = ΔV × (1 − e^{−τ_pulse/τ_m})
    When τ_pulse ≪ τ_m: ≈ ΔV × τ_pulse/τ_m."""
    return dV_mV * (1.0 - np.exp(-tau_pulse_s / tau_m_s))

# ============================================================================
# (A) Single-AP instantaneous peak — crack quasi-steady state (with lateral escape) [the core of this section]
# ============================================================================
# Physics: K⁺ is released at flux rate j_flux = j/τ_release [mol/cm²/s] into a crack of height g, radius R_app
#   disk, with the crack laterally open to the bulk ECS. K⁺ diffuses two-dimensionally within the crack; the lateral escape time
#   τ_lat = R_app²/(4D_free) ≈ 0.033 ms, far shorter than the release duration τ_release ≈ 2 ms
#   (release 60× slower than escape), so the crack is approximately quasi-steady: release rate = lateral diffusion outflow rate.
# Steady-state two-dimensional diffusion equation (with uniform volume source j_flux/g):
#   0 = D·∇²u + j_flux/g , boundary u(R_app) = 0 (bulk)
#   Solution (axisymmetric parabola): u(r) = (j_flux/g)·(R_app² − r²)/(4D)
#   Center peak u_0 = (j_flux/g)·R_app²/(4D_free) = j_flux·R_app²/(4·D_free·g)
# Units: j_flux[mol/cm²/s]·R_app²[cm²]/(D[cm²/s]·g[cm]) = [mol/cm³]
print("\n" + "=" * 84)
print("(A) Single-AP instantaneous peak — crack quasi-steady state (with lateral escape) [the core correction here]")
print("    u_peak = j_flux · R_app² / (4 · D_free · g),   j_flux = j / τ_release")
print("=" * 84)
j_flux = j / tau_pulse              # release-period flux rate [mol/cm²/s] (τ_release = 2 ms)
u_A = j_flux * R_app**2 / (4.0 * D_free * gaps_cm)   # [mol/cm³]
print(f"\n  K⁺ efflux rate j_flux = {j_flux:.3e} mol/cm²/s  (j/τ_release, τ_release={tau_pulse*1e3:.0f} ms)")
print(f"  Lateral escape time τ_lat = R_app²/(4D) = {R_app**2/(4*D_free)*1e6:.2f} µs  (≪ τ_release 2 ms → quasi-steady holds)")
print(f"  {'gap [nm]':>8} | {'u_peak [µM]':>12} | {'ΔE_K [mV]':>10} | {'ΔE_K filtered [mV]':>15}")
for g_nm, u in zip(gaps_nm, u_A):
    u_uM = u * MOL_CM3_TO_uM            # mol/cm³ → µM
    dV   = nernst_mV(u * MOL_CM3_TO_mM) # mol/cm³ → mM
    dV_f = integrator_filter(dV, tau_pulse, tau_m)
    print(f"  {g_nm:8.0f} | {u_uM:12.2f} | {dV:10.4f} | {dV_f:15.4f}")
u_A_50 = u_A[1]                          # g=50 nm reference value [mol/cm³]
print(f"\n  → instantaneous peak ~1–9 µM (decreasing as gap increases), Nernst depolarization ~0.01–0.08 mV,")
print(f"    after filtering (τ_m=110 ms integration) ~0.0002–0.0014 mV.")

# Release-time-course sensitivity: τ_release taken as 0.5–2 ms (repolarization K⁺ current front-loaded, peak current concentrated in
#   the early repolarization phase), instantaneous peak ∝ 1/τ_release, correspondingly amplified 1–4×.
print(f"\n  Release-time-course sensitivity (g=50 nm, τ_release taken as 0.5/1.0/2.0 ms):")
for tr in [0.5e-3, 1.0e-3, 2.0e-3]:
    u_tr = (j / tr) * R_app**2 / (4.0 * D_free * gaps_cm[1]) * MOL_CM3_TO_uM
    print(f"    τ_release={tr*1e3:4.1f} ms → u_peak(g=50nm) ≈ {u_tr:5.1f} µM")
print("  → reasonable instantaneous-peak range ~1–34 µM (τ_release 0.5–2 ms, gap 20–150 nm), always at the µM level.")

# ============================================================================
# (B) Single-AP instantaneous peak upper bound — instantaneous constraint (no escape) [non-physical upper bound]
# ============================================================================
# If the total single-AP efflux j [mol/cm²] were all instantaneously confined within a crack of height g (assuming zero escape,
# i.e., release is a δ function and measurement precedes lateral escape), then:
#     u_inst = j / g   [mol/cm² / cm = mol/cm³]
# Unit check: j[mol/cm²]/g[cm] = mol/cm³. This value is reachable only when τ_release ≪ τ_lat;
# in reality τ_release(2 ms) ≫ τ_lat(0.03 ms), so it is a **non-physical upper bound**, used only to anchor the upper bound.
print("\n" + "=" * 84)
print("(B) Single-AP instantaneous peak upper bound — instantaneous constraint (no lateral escape, non-physical upper bound)")
print("    u_inst = j / g")
print("=" * 84)
u_B = j / gaps_cm                 # [mol/cm³]
print(f"\n  {'gap [nm]':>8} | {'u_inst [mM]':>12} | {'ΔE_K [mV]':>10}")
for g_nm, u in zip(gaps_nm, u_B):
    u_mM = u * MOL_CM3_TO_mM       # mol/cm³ → mM
    dV   = nernst_mV(u_mM)
    print(f"  {g_nm:8.0f} | {u_mM:12.3f} | {dV:10.2f}")
print(f"\n  → upper bound ~0.07–0.52 mM (ΔE_K ~0.6–4.3 mV), still below ΔV_th = {dV_th:.1f} mV.")
print("    this value holds only when the crack is fully sealed (no lateral escape); in reality the crack is open,")
print("    escape (0.03 ms) is far faster than release (2 ms), so (B) overestimates by ~τ_release/τ_lat ≈ 60× and is unusable.")

# ============================================================================
# (C) 1 Hz continuous-discharge bulk steady state — cylindrical line source (S6's old value, mislabeled as "near field")
# ============================================================================
# Formula: steady-state solution of the cylindrically symmetric diffusion equation with first-order clearance (Frankenhaeuser & Hodgkin 1956,
#   J Physiol 131:341, PMID 13320339 [c89] diffusion model; Nicholson series reviews):
#     u(r) = (S / 2πD) · K0(r / λ)，  λ = sqrt(D / k)
# where the line-source strength S = Q_b · f / d_s [mol/cm/s] (K⁺ source per cm of axon length per unit time),
#   f = firing frequency [1/s], d_s = varicosity spacing [cm].
# This is the time average (steady state) of the extra-bundle bulk mean field under 1 Hz continuous discharge, **not** the single-AP instantaneous peak.
print("\n" + "=" * 84)
print("(C) 1 Hz continuous-discharge bulk steady state — cylindrical line source (= S6's old 'near field' value)")
print("    u(r) = (S/2πD)·K0(r/λ),  S = Q_b·f/d_s,  λ = sqrt(D/k)")
print("=" * 84)
f    = 1.0                        # basal firing frequency [1/s]
d_s  = 5e-4                       # varicosity spacing [cm] = 5 µm (mid-value of Goyal & Chaudhury 2013 2–10 µm)
S_line = Q_b * f / d_s            # line-source strength [mol/cm/s]
print(f"\n  Total K⁺ efflux per varicosity Q_b = {Q_b:.3e} mol/AP")
print(f"  Line-source strength S = Q_b·f/d_s = {S_line:.3e} mol/cm/s")
for kname, kk in [("Förstl 1982 (k=0.04/s)", k_fast), ("Hoppe 1991 (k=0.007/s)", k_slow)]:
    for Dname, D in [("D_free", D_free), ("D_tort", D_tort)]:
        lam = np.sqrt(D / kk)
        r   = b_v                      # varicosity surface (outside the terminal readout point ~20–150 nm, same r≈a)
        u   = (S_line / (2 * np.pi * D)) * k0(r / lam)   # [mol/cm³]
        u_nM = u * MOL_CM3_TO_nM       # mol/cm³ → nM
        dV   = nernst_mV(u * MOL_CM3_TO_mM)
        print(f"  {kname:22s} | {Dname:7s} | λ={lam*1e4:5.1f} µm | "
              f"u(r=a)={u_nM:6.2f} nM | ΔE_K={dV:.2e} mV")
print(f"\n  → bulk steady state ~6–7 nM (ΔE_K ~5e-5 mV), about 3 orders of magnitude smaller than the instantaneous peak (A).")
print("    this is the source of S6's original 'terminal near field ~1–7 nM' — actually a bulk time average, not the terminal near-field instantaneous value.")

# ============================================================================
# Summary conclusions (citing the actually computed values, not hard-coded numbers)
# ============================================================================
u_A_uM   = u_A_50 * MOL_CM3_TO_uM       # g=50 nm instantaneous peak [µM]
dA_raw   = nernst_mV(u_A_50 * MOL_CM3_TO_mM)          # g=50 nm instantaneous peak raw depolarization [mV]
dA_filt  = integrator_filter(dA_raw, tau_pulse, tau_m) # after filtering [mV]
gap_raw  = dV_th / dA_raw               # raw gap factor (multiplicative)
gap_filt = dV_th / dA_filt              # gap factor after filtering
# Δ[K⁺]o required for Nav1.7 5 mV subthreshold amplification (Nernst inverse: dV = NERNST·ln((K0+x)/K0) → x)
dK_Nav17 = K_rest_mM * (np.exp(dV_Nav17 / NERNST) - 1.0)   # [mM] ≈ 0.82 mM

print("\n" + "=" * 84)
print("Summary conclusions")
print("=" * 84)
print(f"""
  ① Is a steady state required for activation? — No.
     the terminal reads the single-AP instantaneous peak (A) ~{min(u_A)*MOL_CM3_TO_uM:.0f}–{max(u_A)*MOL_CM3_TO_uM:.0f} µM
     (g=20–150 nm), not the bulk steady state (C) ~6 nM. The two differ by ~3 orders of magnitude, but both are far below the firing threshold.
     baseline [K⁺]o = 4 mM, so a near-field instantaneous value of ~µM is a reasonable magnitude (user's intuition confirmed); ~nM is the single-fiber
     1 Hz time average, not the instantaneous value the terminal reads.

  ② Depolarization of the three concentrations (raw / after the τ_m=110 ms integrator filter):
     (A) instantaneous peak g=50 nm {u_A_uM:.1f} µM  → {dA_raw:.3f} mV / {dA_filt:.4f} mV after filtering
     (B) instantaneous upper bound ~0.07–0.52 mM → ~0.6–4.3 mV (non-physical upper bound)
     (C) bulk steady state ~6 nM       → ~5e-5 mV
     threshold voltage ΔV_th = {dV_th:.1f} mV (derived from the empirical u_th=6.8 mM; Type-A reference RMP {V_rest}→threshold {V_thr} mV, [c134])

  ③ Conclusion: the single-varicosity K⁺ flux (whether the instantaneous peak or the steady state) cannot reach the firing threshold.
     - instantaneous peak (A) g=50 nm is ~{gap_raw:.0f}× from the AP threshold (~{dV_th:.1f} mV) (raw); ~{gap_filt:.0f}× after the integrator filter.
     - even using Nav1.7's lower bound (~5 mV subthreshold amplification), it still requires Δ[K⁺]o ≈ {dK_Nav17:.2f} mM, still ~µM instantaneous
       peak is ~2–3 orders of magnitude higher → the flux cannot even reach the "minimal detectable depolarization".
     - for the flux to reach the mM threshold, an "inventory event" is required: many cells releasing synchronously and spatially confined (pathological high resistance /
       needle-source damage), which is the footing of the three effects (dam-break / resonance) in main text §2.3.5 — the single-varicosity flux is never the source.

  ④ Correction to S6: change 'terminal near field ~1–7 nM' to
     'single-AP instantaneous peak ~1–9 µM (crack quasi-steady state with lateral escape); bulk steady state ~6 nM'.
     the conclusion is unchanged, but the near-field magnitude is corrected from nM to µM (instantaneous), and the gap is corrected from '4–6 orders of magnitude' to
     'instantaneous peak ~2–3 orders of magnitude (raw) / ~4–5 orders of magnitude after integration; bulk steady state ~6 orders of magnitude'.
""")

# ============================================================================
# (D) Complete flux (adding source b VSMC + the endothelial non-source geometric negative result — (A)(B)(C) above only computed source-a varicosities)
# ============================================================================
# S6's "three K⁺ sources":
#   (a) sympathetic varicosity AP → K⁺ (repolarizing Kv) — computed above, instantaneous ~1–9 µM
#   (b) vascular smooth muscle (VSMC) relaxation → K⁺ (BKCa) — added in this section
#   endothelial EDH (SK/IK, non-source, already moved out of the source table) — added in this section (Edwards 1998 myoendothelial gap mM)
#   (c) cell damage → K⁺ (only near the needle puncture) — dam-break-triggered, reaches mM, not included in the resting flux sum
print("\n" + "=" * 84)
print("(D) Complete flux (adding source b VSMC + the endothelial non-source geometric negative result)")
print("=" * 84)

# --- Source (b): VSMC relaxation K⁺ efflux (Ca²⁺ spark → BKCa → STOC circulating current) ---
# [2026-08-22 correction, aligned with model_source_b.py] the old model j_VSMC = C_m·ΔV_hyp/F
#   (single 25 mV hyperpolarization net charge = capacitive lower bound 0.26 pmol/cm²) **underestimates by ~22×**: the K⁺ efflux is a BKCa
#   tonic **circulating current** (continuously offsetting the inward leak current, not a one-shot movement), of magnitude = STOC circulating current.
#   ΔV_hyp = 25 mV (BKCa-mediated hyperpolarization amplitude; EDH typically 10–30 mV, Nelson & Quayle 1995,
#   PMID 7733230 [c100]). I_STOC = 250 pA × 100 ms (Benham & Bolton 1986, PMID 2442353
#   [c163] maximal STOC) → per STOC Q_STOC = 25 pC → j_STOC = Q_STOC/(F·S_VSMC) ≈ 5.6 pmol/cm².
dV_hyp  = 25e-3                     # [V] = 25 mV
j_hyp   = C_m * dV_hyp / F          # [mol/cm²] capacitive lower bound (single 25 mV net charge) = 0.26 pmol/cm²
# VSMC surface area (Smirnov & Aaronson 1992, PMID 1284312, human mesenteric artery VSMC: capacitance ~46 pF →
#   46e-12 F / 1e-2 F/m² = 4600 µm²; direct measurement, not a geometric lower bound. The old 1131 µm² underestimated by ~4×):
S_VSMC  = 4.6e-5                    # surface area [cm²] = 4600 µm²
I_STOC  = 250e-12                   # [A] STOC peak current (Benham & Bolton 1986 [c163], maximal)
T_STOC  = 100e-3                    # [s] STOC duration
Q_STOC  = I_STOC * T_STOC           # [C] = 2.5e-11 C = 25 pC
j_STOC  = Q_STOC / (F * S_VSMC)     # [mol/cm²] moved per STOC = 5.6 pmol/cm² ≈ 22× the capacitive lower bound
print(f"\n  Source (b) VSMC: STOC circulating current (Ca²⁺ spark→BKCa), {Q_STOC*1e12:.0f} pC per STOC → "
      f"{j_STOC*1e12:.2f} pmol/cm²")
print(f"         = {j_STOC/j_hyp:.0f}× of the capacitive lower bound C_m·ΔV_hyp/F = {j_hyp*1e12:.2f} pmol/cm²"
      f"(circulating current vs one-shot movement)")
print(f"         VSMC surface area ≈ {S_VSMC*1e8:.0f} µm² (human mesenteric 46 pF), total K⁺ per STOC = "
      f"{Q_STOC/F:.2e} mol = {Q_STOC/F/Q_b:.0f}× the varicosity Q_b")
# Source-b readout distance (2026-08-22 measurement-anchored, Luff 2000 [c167] PMID 10602088): source-b K⁺ is released from
#   the **VSMC surface**, and the sensory terminal = **SP+ sensory C-fiber varicosity**, whose distance from the VSMC = **0.4–2 µm**
#   (the majority; SP- sympathetic varicosities are closer <1 µm; SP+ does not form junctions with the SMC → µm-scale diffusion, not a nm crack),
#   **not the source-a 50 nm crack** (that is the sympathetic varicosity↔sensory terminal nerve-nerve apposition [c49]).
#   The readout point is scanned over 0.4/1/2 µm, representative value 1 µm (mid-value of the measured range).
g_b_smc = 1.0e-4                     # [cm] = 1 µm source-b readout point (VSMC surface→SP+ sensory varicosity, Luff 2000 [c167])
# The concentration at this 1 µm readout point uses the full cylindrical-source reaction-diffusion field (model_source_b.py):
#   STOC peak (J_stoc=56 pmol/cm²/s, uniform cylindrical aperture) = 99 µM; transient point-source near-cluster upper bound (1 µm) = 865 µM.
# the old 'u_VSMC_inst = j_STOC/g(50nm)' mistook source b for a 50 nm apposed crack source, a geometric error (overestimating ~2.6×).

# --- Endothelium (non-source): SK/IK → K⁺ (EDH) ---
# Edwards 1998 (Nature 396:269, PMID 9834033 [c125]): ACh → endothelial SKCa/IKCa → K⁺ efflux,
#   myoendothelial gap [K⁺] 4.6→11.6 mM (**myoendothelial gap local value**, mM level).
# Mechanism review [c123] Garland & Dora 2017 (PMID 26752699): endothelial [Ca²⁺]i↑→SKCa/IKCa→hyperpolarization→K⁺ efflux.
# Key geography: the myoendothelial gap is **subendothelial** (intimal side), the sensory terminal is in the **adventitia/media** ([c49] apposed to varicosities,
#   VSMC), separated by the full thickness of the media (arterioles ~several–tens of µm). So the mM myoendothelial value
#   after media diffusion/clearance dilution, is already at the ~µM level when it reaches the sensory terminal (same order as (b), not mM).
print(f"\n  Endothelium (non-source) EDH: Edwards 1998 myoendothelial gap [K⁺] 4.6→11.6 mM (**myoendothelial local mM**,"
      f"[c125] PMID 9834033)")
print("         but the sensory terminal is in the adventitia/media ([c49]) and does not directly read the myoendothelial gap → by the time it reaches the terminal it has been diluted to ~µM.")

# --- Complete flux sum (a + b; endothelium of the same order after dilution) ---
u_a_mid = u_A[1] * MOL_CM3_TO_uM            # varicosity instantaneous g=50nm [µM] ≈ 3.4
u_b_eff = 99.0                              # VSMC 1µm readout-point STOC peak [µM] (model_source_b.py cylindrical aperture)
#   Geometric correction (2026-08-22, Luff 2000 [c167]): the old u_b_eff = 8.6 × 22 ≈ 187 µM used the g=50nm **crack**
#   (the source-a aperture), which is wrong — source-b K⁺ exits from the VSMC surface, the SP+ sensory varicosity is 0.4–2 µm from the VSMC,
#   so the readout point should be 1 µm (cylindrical aperture 99 µM). Transient point-source near-cluster upper bound (1 µm) ≈ 865 µM, see full_field [3.5].
u_total = u_a_mid + u_b_eff                  # total near-field flux [µM] ≈ 102
# Threshold (single aperture): u_th = 6.8 mM (Tkachenko 2023 empirical concentration anchor, total [K⁺]ₒ 10.8 − K_rest)
u_th = 6.8e3                                 # full threshold 6.8 mM → ΔV_th = 26.71·ln(10.8/4.0) ≈ 26.5 mV [µM]
print(f"\n  Complete near-field flux (a+b: source-a varicosity g=50nm crack, source-b VSMC 1µm readout point [Luff 2000 0.4–2 µm]):"
      f"varicosity ~{u_a_mid:.1f} µM + VSMC ~{u_b_eff:.1f} µM ≈ {u_total:.0f} µM")
print(f"  → about {u_th/u_total:.0f}× from the full threshold 6.8 mM (ΔV_th≈26.5 mV).")
print(f"  (the old C·ΔV aperture 'ordinary skin ~170×' is deprecated; STOC circulation 22× + 1µm readout point"
      f"after correction source b ~{u_b_eff:.0f} µM, plus varicosity ~{u_total:.0f} µM.)")
print(f"  → conclusion (direction unchanged): after adding source (b) (endothelium is a non-source), the complete flux is at the ~{u_total:.0f} µM level,"
      f"still about {u_th/u_total:.0f}× below the full threshold,")
print("     still cannot reach the mM threshold; crossing the threshold still relies on an 'inventory event' (dam-break / resonance).")

# ============================================================================
# (E) Along-meridian density: near-field scaling (user's correction — density = varicosity density↑ + mean distance↓, scaling the near field, not bulk)
# ============================================================================
# User's correction: why does the density increase? It stems from increased sympathetic + vascular supply capacity, concretely manifest as
#   ① varicosity density increase → the number of varicosities each sensory terminal faces, N_boutons ↑ (near-field source count ↑)
#   ② mean distance decrease → but **the two distances are separated** (user's 2026-08-22 re-correction):
#       source a = varicosity↔sensory terminal crack g_a (20–150 nm [c49]); source b = VSMC↔sensory terminal readout point g_b
#       (0.4–2 µm, Luff 2000 [c167] measurement) — both are near-field ∝ 1/g, but g_a ≠ g_b,
#       so the scaling cannot share the same 'gap g'.
# Both are **near-field** effects (the terminal reads the sources in its own crack), not (just) a bulk regional average.
# The previous version treated 'density' as bulk scaling, which was wrong — density first changes the terminal's **local coupling** (source count + distance),
#   not the regional total. Near field ∝ N_boutons × (1/g), both rising with density.
print("\n" + "=" * 84)
print("(E) Along-meridian density: near-field scaling (density = varicosity count↑ + gap↓, not bulk)")
print("=" * 84)

# --- Anchored factors (with literature) ---
# Density (2026-08-23 correction, the old 4.8× is deprecated):
#   old 4.8× = acupoint PCV occurrence rate 87.5% / non-acupoint 18.1% (Álvarez-Prats 2019, PMID 31015854),
#   but occurrence rate ≠ density (Poisson ≈10×, clustering <4.8, direction unknown), and only acupoints vs 4cm control points were measured (not along-meridian segments).
#   Split into two coefficients:
#   · along-meridian 1.5× = along-meridian average (acupoint 3.0 × acupoint along-meridian fraction ~0.15–0.33 + inter-acupoint segment ~1.0)
#   · acupoint 3.0× = acupoint local microvascular density (conservative placeholder; Liu 2014/Hsiu 2007 only "significantly higher" with no factor;
#     no difference among Ren meridian CV3/4/5/7 PMC8099720 → no additional amplification in along-meridian segments)
#   Also: vascularization proxy; varicosity ∝ vessel is still an assumption (user: density = sympathetic + vascular supply capacity↑ → varicosity density↑).
density_x  = 1.5    # along-meridian coefficient (main value, PSC conduction along the meridian)
dens_point = 3.0    # acupoint coefficient (local upper bound at the needling point, output for comparison only)
# Recruitment: attention/intention recruits additional varicosities, lower bound ~1.5× (PMID 31347763), VA 0→90% (PMID 15049717).
recruit_x = 2.0

# --- Near-field baseline (source-a varicosity + source-b VSMC, ordinary skin) ---
# Two distances (user's 2026-08-22 correction): source a = varicosity↔sensory terminal crack 20–150 nm [c49];
#   source b = VSMC surface↔SP+ sensory varicosity 0.4–2 µm (Luff 2000 [c167] measurement, representative 1 µm)
#   — **not the same g**.
u_a_mid   = u_A[1] * MOL_CM3_TO_uM   # single varicosity g=50nm ~3.4 µM (Section A)
u_b_eff   = 99.0                      # single VSMC 1µm readout-point STOC peak (cylindrical aperture, see section (D) / full_field)
u_near_base = u_a_mid + u_b_eff       # ≈ 102 µM (the old g=50nm aperture 190 µM is deprecated)

# --- Near-field scaling (density → more varicosities/vessels; recruitment → more varicosity firing) ---
# Source (a) varicosity: × density (more varicosities facing the terminal) × recruitment (more varicosities firing)
u_a_dens = u_a_mid * density_x * recruit_x
# Source (b) VSMC: × density (more vessels); recruitment is a varicosity-firing event and does not directly multiply the VSMC (conservative, no downstream NA amplification added)
u_b_dens = u_b_eff * density_x
u_near_dens = u_a_dens + u_b_dens     # near field (density + recruitment)
u_near_point = u_a_mid * dens_point * recruit_x + u_b_eff * dens_point   # acupoint 3.0× local upper bound
print(f"  Ordinary-skin near field (a+b: source-a g=50nm crack + source-b 1µm readout point [Luff 2000 0.4–2 µm]):"
      f"{u_near_base:.0f} µM (varicosity {u_a_mid:.1f} + VSMC {u_b_eff:.1f})")
print(f"  density {density_x:.1f}× (varicosities+vessels) + recruitment {recruit_x:.1f}× (varicosity firing) → "
      f"near field ≈ {u_near_dens:.0f} µM")
print(f"  acupoint {dens_point:.1f}× (local upper bound): near field ≈ {u_near_point:.0f} µM (needling point, not the along-meridian conduction main value)")
print(f"  → about {u_th/u_near_dens:.1f}× from the full threshold 6.8 mM (near field, raw, anchored part).")

# --- Unanchored factors (qualitative): distance shortening → near field ∝ 1/g ↑ (**the two distances separated**) ---
# Source-a crack g_a: varicosity↔sensory terminal 150→20 nm = 7.5× (20 nm = minimal synaptic cleft, physical upper bound)
# Source-b readout point g_b: VSMC↔SP+ sensory varicosity baseline 0.4–2 µm (Luff 2000 [c167] already anchored); but
#   'is the meridian vs non-acupoint VSMC→sensory terminal distance shorter' still lacks a single PMID (needs meridian-specific site EM quantification).
# TODO [to check]: anatomical quantification of the acupoint/meridian 'VSMC→sensory terminal mean distance' vs non-acupoint, to anchor the g_b along-meridian scaling.
gap_a_max     = 150.0 / 20.0          # 7.5× multiplies only source a (crack); source-b along-meridian shortening is unanchored, not multiplied for now
u_near_gapmax = u_a_dens * gap_a_max + u_b_dens
print(f"\n  Unanchored: source-a crack g_a 150→20 nm = {gap_a_max:.1f}× (multiplies only the source-a varicosity);"
      f"source-b readout point g_b baseline is already anchored (Luff 2000 0.4–2 µm), along-meridian shortening unanchored, not multiplied for now.")
print(f"    if included (most favorable), near field ≈ {u_near_gapmax:.0f} µM → still ~{u_th/u_near_gapmax:.1f}× from the full threshold 6.8 mM.")
print(f"  → anchored conclusion: density×recruitment reduces the near-field gap from ~{u_th/u_near_base:.0f}× (ordinary skin) to ~{u_th/u_near_dens:.1f}× (along-meridian 1.5×, source b already dominant).")
print(f"    source-a crack narrowing (unanchored, upper bound {gap_a_max:.1f}×) reduces it at most to ~{u_th/u_near_gapmax:.2f}×; the four apertures")
print(f"    (ordinary {u_near_base:.0f} / along-meridian {u_near_dens:.0f} / acupoint {u_near_point:.0f} / crack upper bound {u_near_gapmax:.0f} µM) all fail to cross the full threshold 6.8 mM; still need a literature anchor for 'source-a crack along-meridian shortening' (the source-b readout point is already anchored by Luff 2000).")
print(f"  → conclusion (closer to threshold after the STOC circulation correction, direction unchanged): density×recruitment is a **near-field** scaling (source count + distance),")
print(f"    source b (~{u_b_dens:.0f} µM) dominates, along-meridian ~{u_near_dens:.0f} µM is still ~{u_th/u_near_dens:.1f}× from the full threshold 6.8 mM;")
print(f"    all near-field apertures fail to cross the threshold; crossing the full threshold still relies on an 'inventory event' (dam-break / resonance).")
# [c151] Dodge 1994 (PMID 7923638) partially anchors the 'sympathetic→vessel(SMC)' leg of gap_factor:
#   rabbit ear artery (peripheral) vs middle cerebral artery (central) EM quantification: releasable varicosities (bare varicosity, <1 µm from SMC)
#   density = 18 vs 2.8 per 500 µm circumference (peripheral 6.4×); varicosity-SMC mean gap = 1.18±0.06 vs 4.95±0.23 µm
#   (peripheral 4.2× tighter); <500 nm close contact 25% vs 3% (peripheral 8×). = peripheral neurovascular coupling 'denser + closer'.
#   ⚠️ only anchors the 'sympathetic→vessel(SMC)' leg (arterial NEJ), not the 'sensory→sympathetic' leg (sensory terminal↔sympathetic varicosity gap,
#   still no single PMID); and it is the rabbit ear artery (not human skin meridian / not capillary). But it gives a quantitative precedent for 'tighter peripherally'.
#   Note the scale: the varicosity-SMC gap 1.18 µm is one order of magnitude larger than the sensory terminal crack 20–150 nm (NA diffusion µm-scale vs synaptic K⁺ nm-scale).
print(f"  [c151] Dodge 1994 (PMID 7923638) partially anchors the gap_factor 'sympathetic→vessel' leg: releasable varicosity density")
print(f"    18 vs 2.8 per 500 µm circumference (peripheral 6.4×), varicosity-SMC gap 1.18 vs 4.95 µm (peripheral 4.2× tighter),")
print(f"    <500 nm close contact 25% vs 3% (peripheral 8×) = EM quantitative precedent for peripheral neurovascular coupling 'denser + closer'.")
print(f"    ⚠️ the 'sensory→sympathetic' leg is still missing (sensory terminal↔sympathetic varicosity gap, no single PMID); and it is the rabbit ear artery, not the skin meridian.")

# --- Bulk dimension of the density increase (user's addition): release capacity↑ + clearance↓ = precondition for an inventory event ---
# The previous section treated density as near-field scaling (varicosity count↑ + gap↓). User's addition: density also has a **bulk** dimension —
#   density↑ = ① local total K⁺ release capacity↑ (source count per volume ∝ vessel/varicosity density)
#          ② more cells/vessels/structures → ECS tortuosity λ↑ → effective diffusion D_eff=D/λ² ↓
#             → diffusion clearance pressure↓ (K⁺ harder to escape, retention↑).
#   This is exactly the **precondition** for an 'inventory event' (dam-break / resonance): turning the local tissue into a better K⁺ reservoir (more sources + slower leak).
#   The difference from dam-break is magnitude: the density-induced λ rise is small (unanchored, estimated ~1.3–1.6×),
#   while the dam-break clearance shutdown is large (k 4.6→0.4×10⁻³/s, ~1000×, Syková ischemia).
# Magnitude: density 1.5× (along-meridian) / 3.0× (acupoint) × clearance ~1.5–2× (λ↑) ≈ 2–6× bulk scaling → does not close the µM→mM gap (~5 orders);
#   closing still relies on dam-break (clearance shutdown) or resonance (whole-cell K⁺ store synchronization). Density is the precondition that 'lowers the gate', not the 'gate closure' itself.
print(f"\n  [user's addition · bulk dimension] density↑ also has a bulk effect: ① release capacity↑ (source density ∝ vessel/varicosity)")
print(f"    ② clearance↓ (more structures→λ↑→D_eff↓). = inventory-event precondition (more sources + slower leak), same physics as dam-break but smaller magnitude:")
print(f"    density ~1.5× (along-meridian) / ~3.0× (acupoint) × clearance ~1.5–2× ≈ 2–6× bulk scaling, does not close the µM→mM gap (~5 orders).")
print(f"    ⚠️ the density-induced λ rise is unanchored (no direct measurement of acupoint tortuosity) — this is an estimate, to be supplemented.")

# ============================================================================
# (F) Capillary wall: K⁺ sensing/conduction, not a (main) source (user's follow-up: how much K⁺ does the capillary wall supply during the vasodilation phase)
# ============================================================================
# Literature (verified online, three complementary items):
#   (i)   Jackson 2017 Microcirculation 24:e12319, PMID 27652592 — authoritative review of endothelial Kir,
#         providing three quantitative anchors (A1–A3, see below).
#   (ii)  Longden 2017 Nat Neurosci 20:717, PMID 28319610 — active-neuron AP repolarization releases
#         K⁺ ([K⁺]o → ~10 mM) → sensed by **capillary endothelial Kir2.1** → hyperpolarization → along the endothelium
#         gap junctions (electrical syncytium) **retrograde conduction mm** → upstream arteriole dilation (functional hyperemia).
#   (iii) Longden, Hill-Eubanks & Nelson 2016, J Cereb Blood Flow Metab 36:492,
#         PMID 26661232 — capillary endothelium **lacks SK/IK (Ca²⁺-activated K⁺ channels)**, so it cannot, like
#         the arteriolar endothelium (non-source, Edwards 1998), release K⁺ via SK/IK.
#
# Jackson 2017 three quantitative anchors (PMID 27652592):
#   A1. capillary/microvascular endothelial RMP ≈ −30~−40 mV (main text 'intact pressurized vessels
#       ~30–40 mV') → anchors taking ΔV_cap = 5–10 mV below as reasonable.
#   A2. **only ~144 functional Kir channels per endothelial cell** (mouse mesenteric artery, Sonkusare 2016 data)
#       — very few channels; even if all open, the K⁺ efflux is tiny → the capillary is not an effective K⁺ source.
#   A3. **K⁺ sensing threshold ~8–10 mM** (active-tissue interstitial K⁺ accumulation 'on the order of 8–10 mM';
#       3–5 → 8–15 mM converted to hyperpolarization) — this is the [K⁺]o required to activate the capillary Kir,
#       which is ~1–2 orders of magnitude higher than this model's near-field instantaneous peak (~81 µM ordinary / ~127 µM along-meridian 1.5×).
#       ⇒ dual conclusion: (a) the capillary is a **downstream sensing element** that 'reads mM', not a source that 'produces µM';
#       (b) the near-field µM transient cannot even activate the capillary Kir, corroborating that 'crossing the threshold requires bulk dam-break / resonance'.
#   ⚠️ subtype correction: capillary endothelium expresses Kir2.1 (essential subunit, knockout abolishes activity) + Kir2.2/2.3 mRNA
#      (Jackson 2017 main text); precise wording = 'has Kir2.x, **no SK/IK**', not 'only Kir2.1'.
#
# Quantification (charge conservation, same as source b/c):
#   ΔV_cap = 5–10 mV hyperpolarization (caused by the Kir outward current when [K⁺]o goes 3→15 mM, small amplitude):
#   j_cap = C_m·ΔV_cap/F — the K⁺ efflux required to hyperpolarize the EC by 5–10 mV.
print("\n" + "=" * 84)
print("(F) Capillary wall: K⁺ sensing/conduction, not a (main) source")
print("=" * 84)

dV_cap = np.array([5e-3, 10e-3])        # capillary EC hyperpolarization 5–10 mV (Jackson 2017/Longden 2017, small amplitude)
j_cap  = C_m * dV_cap / F               # K⁺ efflux [mol/cm²] (charge conservation)
print(f"  capillary EC hyperpolarization ΔV_cap = 5–10 mV (RMP −30~−40 mV → more negative; Jackson 2017)")
print(f"  j_cap = C_m·ΔV_cap/F = {j_cap[0]*1e12:.2f}–{j_cap[1]*1e12:.2f} pmol/cm²")
print(f"  → 1/{j/j_cap[0]:.0f}–1/{j/j_cap[1]:.0f} of the varicosity {j*1e12:.2f} pmol/cm²;"
      f"1/{j_STOC/j_cap[0]:.0f}–1/{j_STOC/j_cap[1]:.0f} of the VSMC {j_STOC*1e12:.2f} pmol/cm².")
print(f"  channel-count anchor (Jackson 2017): microvascular EC has only ~144 functional Kir channels/cell,")
print(f"    and the EC high membrane resistance 2–100 GΩ — a few channels suffice to move the membrane potential, showing that the role of Kir is")
print(f"    'membrane-potential modulation' (hyperpolarization amplifier), not massive K⁺ efflux; the reverse proves the capillary is not a K⁺ source.")
# K⁺ sensing-threshold consistency (Jackson 2017 A3): the capillary Kir needs ~8–10 mM, vs near-field µM
K_sense_lo, K_sense_hi = 8.0, 10.0      # [mM] active-tissue interstitial K⁺ accumulation (Jackson 2017 main text)
print(f"  K⁺ sensing threshold ~{K_sense_lo:.0f}–{K_sense_hi:.0f} mM (Jackson 2017) vs near-field instantaneous "
      f"{u_total:.0f} µM (ordinary) / {u_near_dens:.0f} µM (along-meridian):")
print(f"    higher by ~{K_sense_lo*1e3/u_total:.0f}–{K_sense_hi*1e3/u_total:.0f}× (ordinary),"
      f"~{K_sense_lo*1e3/u_near_dens:.0f}–{K_sense_hi*1e3/u_near_dens:.0f}× (along-meridian).")
print(f"    ⇒ the capillary reads the mM level (bulk dam-break domain); the near-field µM cannot even reach the capillary Kir.")
# How is this 8–10 mM reached? (user's follow-up) — bulk synchronous discharge + ceiling, not the single-cell near field
#   Heinemann & Lux 1977, Brain Res 120:231, PMID 832122 [c148]: cat cortex synchronous/repetitive stimulation →
#   extracellular K⁺ rises to the 'ceiling' ≈10.2 mequiv/l (~10 mM), maintained by glial active uptake + reduced release
#   at high K⁺ (driving force↓) reaching a 'release = uptake' steady-state ceiling.
#   Syková 1983 [c84] (PMID 6139844): cortex resting 2.8–3.4 → max 10.2 mM → >12 mM triggers SD.
#   mechanism = **bulk synchronous release** (many cells simultaneously releasing K⁺ during AP repolarization) accumulating faster than clearance
#   (glial Na/K pump + diffusion) in the confined ECS → i.e., this model's 'inventory event' (dam-break / resonance), not the single-varicosity near-field µM.
print(f"  how 8–10 mM is reached: bulk synchronous discharge + ceiling, not the single-cell near field —")
print(f"    Heinemann & Lux 1977 [c148] (PMID 832122): cat cortex synchronous stimulation → K⁺ rises to ceiling")
print(f"    ≈10.2 mequiv/l, balanced by glial active uptake + reduced high-K⁺ release; Syková 1983 [c84]: cortex")
print(f"    2.8–3.4→10.2 mM, >12 mM triggers SD. = this model's 'inventory event' (bulk synchrony + confinement + clearance lag).")
# Key caveat: the ceiling is a CNS (cortex) value; peripheral nerve at 36°C rarely exceeds 1 mM (Hoppe 1991 [c112])
print(f"  ⚠️ but 10–12 mM is a CNS cortex value: peripheral nerve at 36°C rarely exceeds 1 mM (Hoppe 1991 [c112],")
print(f"    temperature-dependent, several mM only at 22°C). The meridian is **peripheral** tissue → the capillary Kir 8–10 mM is even")
print(f"    harder to reach along the meridian (peripheral ceiling ~1 mM, still below the full threshold 6.8 mM) → again confirms 'crossing the threshold requires bulk inventory'.")
# Graded K⁺ vasodilation (user's follow-up: is the Laogong warming a capillary dilation signal):
#   K⁺ vasodilation is **graded**, not only the single 8–10 mM capillary Kir path —
#   Burns & Cohen 2004 [c149] (PMID 15280082): hamster cremaster arteriole K⁺ 5→20 mM dilation,
#   ouabain (Na/K-ATPase block) essentially eliminates the hyperpolarization + dilation, Ba²⁺ (Kir block) only shortens the time course →
#   **Na/K-ATPase is the main driver (low threshold, Km~1–2 mM), Kir is a high-threshold secondary path (time course/kinetics)**.
#   Implication: the peripheral Laogong warming (intention→sympathetic→local K⁺ ~1–2 mM) can be driven by Na/K-ATPase,
#   not requiring the 8–10 mM capillary Kir — 'capillary Kir is unreachable peripherally' ≠ 'no K⁺ vasodilation peripherally'.
print(f"  graded K⁺ vasodilation (Burns & Cohen 2004 [c149], PMID 15280082):")
print(f"    the main K⁺ dilation path is Na/K-ATPase (abolished by ouabain, low threshold Km~1–2 mM); Kir (8–10 mM)")
print(f"    only tunes the time course/kinetics (Ba²⁺ block only shortens the time course). So the peripheral Laogong warming (~1–2 mM) can be driven by Na/K-ATPase")
print(f"    without the 8–10 mM capillary Kir — 'capillary Kir unreachable' ≠ 'no K⁺ vasodilation peripherally'.")
# The warming phenomenon itself has literature, do not mark it 'to check': intention/imagery→skin sympathetic→skin temperature self-regulation =
#   [c72] Taub & Emurian 1976 (PMID 990346): single feedback site, learned in 4 training sessions, mean 2.2°F (up to 6.5°F),
#   9–14°F bidirectional regulation within 15 min after training = intentional-trainability modulation of skin sympathetic vasomotor output (phenomenon anchored).
#   Only the **mechanistic attribution** remains to check: the local [K⁺]o during intention was not directly measured → K⁺-driven vs neurogenic (autonomic) dilation is undifferentiated.
print(f"  the warming phenomenon itself has literature (not 'to check'): [c72] Taub & Emurian 1976 (PMID 990346)")
print(f"    intention/biofeedback→skin temperature self-regulation 2.2°F (up to 6.5°F), 9–14°F/15min after training = intentional-trainability modulation")
print(f"    of skin sympathetic vasomotor output. ⚠️ only the mechanistic attribution remains to check: [K⁺]o was not directly measured during intention →")
print(f"    K⁺-driven (Na/K-ATPase) vs neurogenic dilation (autonomic) is undifferentiated; ouabain/Ba²⁺/autonomic block")
print(f"    or a skin K⁺ microelectrode (sweat contamination must be controlled) is needed for a definitive conclusion.")
# Redundancy caveat (Joyner & Wilkins 2007 [c152], PMID 17640934): vasodilation has **no single obligatory factor** —
#   exercise hyperemia (skeletal-muscle blood flow↑100×) has been studied 100+ years without finding an obligatory factor; NO/adenosine/ATP/KATP single, double, or triple blockade leaves hyperemia intact.
#   K⁺ (Na/K-ATPase/Kir) is one of many redundant factors, not obligatory. Two constraints: ① sympathetic withdrawal is only ~3× (upper bound),
#   ② humans lack cholinergic sympathetic vasodilator nerves. So 'Laogong warming = K⁺-driven' should be stated as 'one candidate mechanism', not 'main path / obligatory'.
#   ⚠️ the redundancy argument instead corroborates this model's multi-source architecture of 'three K⁺ sources + ATP/P2X non-ionic' (multi-source is necessary, not a defect).
print(f"  [c152] Joyner & Wilkins 2007 (PMID 17640934) redundancy caveat: vasodilation has **no single obligatory factor** —")
print(f"    exercise hyperemia has been studied 100+ years without finding an obligatory factor (NO/adenosine/ATP/KATP single, double, or triple blockade leaves hyperemia intact). So K⁺")
print(f"    (Na/K-ATPase/Kir) is **one of** the redundant factors, not obligatory; and sympathetic withdrawal is only ~3×, humans lack")
print(f"    cholinergic sympathetic vasodilation → 'Laogong warming = K⁺-driven' should be stated as 'one candidate mechanism', not 'main path / obligatory'.")
print(f"    ⚠️ the redundancy argument instead corroborates this model's 'three K⁺ sources + ATP/P2X' multi-source architecture — multi-source is necessary, not a defect.")
print(f"  conclusion: the capillary-wall K⁺ efflux ~0.05–0.1 pmol/cm² (only ~144 Kir channels/cell) is")
print(f"    ~1/10–1/20 of the varicosity, and is a **response (amplifier)** to neuronal K⁺ (the true source, ~10 mM), not an independent source.")
print(f"    capillary = K⁺ sensing/electrical conduction element ('K⁺ initiation + electrical conduction', Longden 2017),")
print(f"    **does not constitute a new extracellular K⁺ accumulation source** — not added to the three-source sum.")
# [c150] Poli 2026 (PMID 41609600, cardiac NVU) addition: sympathetic varicosities **directly appose the capillary endothelium** (EM Fig 1C,
#   endothelium expresses COMT) = direct anatomical evidence of the last hop of 'sensory nerve→sympathetic→vessel(capillary)' (cardiac, not skin).
#   But the capillary-side communication is **NA→endothelial/pericyte adrenergic** (synapse-like or paracrine spillover, not distinguished by the authors),
#   not K⁺ — complementary to 'capillary has no SK/IK = non-K⁺ source': the capillary is a sympathetic target, and the terminal effector is adrenergic.
#   ⚠️ the causal direction is counterintuitive: sympathetic **loss**→capillary density↑ + caliber↓ (compensatory remodeling 155.8→196.4/area),
#   so 'sympathetic-vessel density↑' should be stated as 'sympathetic nerve **innervation density**↑', not mistakenly as 'sympathetic→capillary **count**↑'.
print(f"  [c150] Poli 2026 (PMID 41609600, cardiac NVU) addition: sympathetic varicosities EM-directly appose the capillary endothelium")
print(f"    (Fig 1C, endothelium expresses COMT) = direct anatomical evidence of the last hop of 'sensory nerve→sympathetic→vessel(capillary)'.")
print(f"    but the capillary-side communication is NA→endothelial/pericyte adrenergic (not K⁺), complementary to 'capillary has no SK/IK = non-K⁺ source'")
print(f"    : the capillary is a sympathetic target and the terminal effector is adrenergic. ⚠️ the causal direction is counterintuitive — sympathetic loss→capillary")
print(f"    density↑ + caliber↓ (compensatory remodeling 155.8→196.4/area), so 'sympathetic-vessel density↑' should be stated as")
print(f"    'sympathetic nerve innervation density↑', not mistakenly as 'sympathetic→capillary count↑'.")

# ============================================================================
# Appendix: scaling of a whole bundle of ~40 fibers (only sums the bulk, does not change the instantaneous peak read by a single terminal)
# ============================================================================
print("=" * 84)
print("Appendix: scaling of a whole bundle of ~40 fibers firing synchronously")
print("=" * 84)
N_fib = 40
u_bulk_40 = S_line * N_fib / (2 * np.pi * D_free) * k0(b_v / np.sqrt(D_free / k_fast)) * MOL_CM3_TO_nM
print(f"  bulk steady state × 40 fibers ≈ {u_bulk_40:.0f} nM (still nM level, ~6 orders from threshold).")
print("  the instantaneous peak (A) does not scale with fiber count: each sensory terminal reads the crack of the varicosity it faces,")
print("  not 40 varicosities sharing one crack (varicosity spacing 2–10 µm, each independent).")
