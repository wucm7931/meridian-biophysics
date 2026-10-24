# -*- coding: utf-8 -*-
"""
Full reaction-diffusion field recomputation for source b — three-tier ladder closed loop
============================================================

Purpose (after user correction 2026-08-21): close the loop of the full-field recomputation of the "three-tier ladder (terminal integration + clearance + density)" for
**source b** (K⁺ efflux from the vascular smooth muscle BKCa relaxation state), answering:
"Does density × recruitment cross threshold (u_th = 6.8 mM full threshold)?"

Old assumption overturned (corrected head-on by this script):
  The old "source b = 8.6 µM" came from the C·ΔV/F capacitance lower bound of 0.26 pmol/cm² (the net charge of one 25 mV hyperpolarization).
  This is a **conceptual error**: K⁺ efflux is a **sustained circulating current** (each STOC must keep countering
  the inward leak current), not a one-shot transport. A single STOC moves ~5.6 pmol/cm² ≈ 22× the capacitance lower bound.

Physical chain (literature-anchored, see each section): Ca²⁺ spark → BKCa → STOC (250 pA × 100 ms × ~1 Hz)
  → K⁺ efflux → hyperpolarization → relaxation. This is **the relaxation mechanism itself, the tonic primary source** (not NA-driven).

Unit iron rule: 1 mol/m³ ≡ 1 mM (consistent with common.py). Hence
  · µM = value × 1e3; · pmol/cm² = mol/m² × 1e8; · µA/cm² = (mol/m²/s·F) × 100.

Geometry: the arteriole is a finite-radius cylinder (r₀≈12.5 µm); use the cylindrical-source reaction-diffusion solution with first-order clearance k
  u(r) = (J·λ/(D_eff·α)) · K₀(r/λ) / K₁(r₀/λ)   (r ≥ r₀, λ=√(D/k))
  Note: for r ≥ r₀ this solution is approximately equal to the line-source K₀ solution when r₀/λ ≪ 1 (K₁(x)≈1/x),
  so the claim "line source overestimates by 1000×" does not hold — the real errors of the old code were ① using C·ΔV instead of the BKCa
  circulating flux (22×, after changing S_CELL to 4600 µm²), and ② not including the geometric treatment of α (1 µm gap vs the cylindrical λ field).

Conclusion preview (details at the end):
  · Source b is the primary source: unit membrane-area flux = 22×/s of the capacitance lower bound, ~5.6× that of varicosity source a.
  · But concentration is set by the "flux ÷ clearance" steady-state balance: ordinary-skin STOC peak (cylindrical uniform treatment) ~99 µM,
    time-average ~10 µM; along-meridian 1.5× (acupoint 3.0×) → STOC peak ~0.15 mM (acupoint ~0.30 mM),
    time-average ~15 µM (all below the 6.8 mM full threshold).
  · Transient point-source treatment (STOC localized, S-independent) near-cluster peak 0.9 mM (1 µm) – 2.2 mM (0.4 µm),
    physically favoring the latter (localized point source, not whole-cell uniform); but all four treatments are below the full threshold 6.8 mM.
    The gap is compressed from the old ~125× to ~45× (along-meridian 1.5× uniform treatment).

Run: python model_source_b.py
"""

import numpy as np
import sys
from scipy.special import k0, k1

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

print("=" * 86)
print("Full reaction-diffusion field recomputation for source b (cylindrical source · BKCa circulating flux · density · integrator)")
print("=" * 86)

# ============================================================================
# 0. Physical constants and tissue parameters (per-item literature provenance)
# ============================================================================
F       = 96485.0        # Faraday constant [C/mol]
R_gas   = 8.314          # gas constant [J/(mol·K)]
T       = 310.0          # 37°C [K]
NERNST  = R_gas * T / F * 1000.0   # 26.71 mV (natural log of RT/F)

D_free  = 1.9e-9         # m²/s free-water D_KCl (=1.9e-5 cm²/s; Syková & Nicholson 2008, PMID 18923183)
lam_T   = 1.6            # tortuosity (Nicholson & Syková 1998, PMID 9610885 [c99])
D_eff   = D_free / lam_T**2   # = 7.42e-10 m²/s = 742 µm²/s effective diffusion
ALPHA   = 0.30           # ECS volume fraction (nearfield_l2 uses 0.30; Nicholson 0.2 as conservative lower bound)
K_CLEAR = 0.01           # /s first-order clearance (muscle-anchored Nielsen 2003/Lo 2004/Syková 1994, range 0.003–0.04; Förstl 1982 SCG 0.04 upper end)
LAMBDA  = np.sqrt(D_eff / K_CLEAR)   # clearance length [m] ≈ 2.72e-4 m = 272 µm

# Sensory-terminal reach threshold (Vasylyev & Waxman 2012 as A-class transduction-layer reference; Bostock 2003 [c135])
TAU_M   = 110e-3         # s membrane-integrator time constant τ_m = 110 ms (Bostock 2003 [c135])
K_REST  = 4.0            # mM resting [K⁺]o (interstitial fluid 3.5–5 mM [c321], take median)
DK_FULL = round(10.8 - K_REST, 1)  # Δ[K⁺]o [mM] needed to reach threshold = 6.8 (the single threshold, Tkachenko 2023 empirical anchor)
DV_TH   = NERNST * np.log((K_REST + DK_FULL) / K_REST)  # 26.5 mV threshold voltage (derived from the empirical u_th)

print("\n[0] Tissue parameters")
print(f"    D_eff={D_eff*1e12:.0f} µm²/s (D_free/λ_T², λ_T={lam_T}); α={ALPHA}; "
      f"k={K_CLEAR}/s → clearance length λ=√(D/k)={LAMBDA*1e6:.0f} µm")
print(f"    Terminal: reach threshold {DV_TH:.1f} mV (Δ[K⁺]o≈{DK_FULL:.1f} mM)")

# ============================================================================
# 1. Source-b flux (BKCa circulating current, not the C·ΔV capacitance lower bound) — the core correction of this script
# ============================================================================
# STOC (Benham & Bolton 1986, PMID 2442353): amplitude ~250 pA, duration ~100 ms,
#   75–100 BKCa channels; frequency 0.5–1 Hz/cell (Jaggar 2000). K⁺ efflux is **circulating**
#   (each STOC counters the inward leak current), not a one-shot transport.
I_STOC      = 250e-12     # A = 250 pA (Benham & Bolton 1986, V_m≈−40 mV)
T_STOC      = 100e-3      # s = 100 ms
F_STOC      = 1.0         # Hz baseline frequency (Jaggar 2000 Table1 unstimulated 0.5–1/s)
S_CELL      = 4600e-12    # m² = 4600 µm² single smooth-muscle-cell surface area (Smirnov & Aaronson 1992, PMID 1284312,
#   human mesenteric-artery VSMC: whole-cell-clamp capacitance ~46 pF → 46e-12 F / 1e-2 F/m² = 4600 µm²; normally distributed,
#   constant across 0.4–3.6 mm vessels). The old 1131 µm² was a provenance-free geometric lower bound, underestimating the measurement by ~4× → flux overestimated ~4×.
#   Cross-species robust range (not entering the main value): McGahon 2006 (PMID 17040965, rat retinal arteriole) capacitance 13.1 pF→1310 /
#   confocal ellipsoid→1220 µm²; O'Dwyer 2020 (PMID 32015129, mouse mesentery) 192,000 CaV1.2÷80/µm²→~2400 µm².
J_STOC_MOL  = I_STOC * T_STOC / (F * S_CELL)   # single-STOC K⁺ efflux [mol/m²] = 5.6 pmol/cm²
J_AVG_MOL   = J_STOC_MOL * F_STOC              # time-average flux [mol/m²/s] = 5.6 pmol/cm²/s
J_STOC_flux = I_STOC / (F * S_CELL)            # instantaneous flux during STOC [mol/m²/s] = 56 pmol/cm²/s
# C·ΔV capacitance lower bound (old value, for comparison):
CM          = 1e-2        # F/m² = 1 µF/cm² specific membrane capacitance
DV_HYP      = 25e-3       # V = 25 mV relaxation hyperpolarization amplitude
Q_CAP_MOL   = CM * DV_HYP / F                 # = 2.59e-9 mol/m² = 0.259 pmol/cm²
RATIO_STOC  = J_STOC_MOL / Q_CAP_MOL          # ≈ 22× (single STOC vs capacitance lower bound)
RATIO_AVG   = J_AVG_MOL / Q_CAP_MOL           # ≈ 22×/s (time-average vs capacitance lower bound)

print("\n[1] Source-b flux (BKCa circulating current, not the C·ΔV capacitance lower bound)")
print(f"    STOC: {I_STOC*1e12:.0f} pA × {T_STOC*1e3:.0f} ms × {F_STOC} Hz = "
      f"{I_STOC*T_STOC*1e12:.0f} pC/STOC")
print(f"    Single-STOC K⁺ efflux = {J_STOC_MOL*1e8:.1f} pmol/cm²"
      f" ({RATIO_STOC:.0f}× the capacitance lower bound {Q_CAP_MOL*1e8:.2f} pmol/cm²)")
print(f"    Time-average flux = {J_AVG_MOL*1e8:.1f} pmol/cm²/s = {J_AVG_MOL*F*100:.1f} µA/cm²"
      f" ({RATIO_AVG:.0f}× the capacitance lower bound per second)")

# ============================================================================
# 2. Cylindrical-source reaction-diffusion steady state (finite radius r₀, with α)
# ============================================================================
R_VESSEL = 12.5e-6       # m terminal arteriole radius (d=25 µm; range 10–25 µm)

def cyl_u(r, J, alpha=ALPHA):
    """Cylindrical-source steady state [mol/m³ ≡ mM]: u(r) = (J·λ/(D·α))·K₀(r/λ)/K₁(r₀/λ), r ≥ r₀."""
    pref = J * LAMBDA / (D_eff * alpha)
    return pref * k0(r / LAMBDA) / k1(R_VESSEL / LAMBDA)

print("\n[2] Cylindrical-source steady state (arteriole r₀ = %.1f µm; r₀/λ = %.4f ≪ 1)"
      % (R_VESSEL * 1e6, R_VESSEL / LAMBDA))
u_tonic_r0 = cyl_u(R_VESSEL, J_AVG_MOL)          # time-average steady-state surface concentration
u_stoc_r0  = cyl_u(R_VESSEL, J_STOC_flux)        # quasi-steady surface concentration during STOC
print(f"    Time-average flux J_avg={J_AVG_MOL*1e8:.1f} pmol/cm²/s → surface steady state "
      f"u(r₀)={u_tonic_r0*1e3:.1f} µM")
print(f"    STOC-period flux J_stoc={J_STOC_flux*1e8:.0f} pmol/cm²/s → surface quasi-steady state "
      f"u(r₀)={u_stoc_r0*1e3:.0f} µM")
print(f"    (Note: the old '8.6 µM' = C·ΔV one-shot transport; this solution = the steady-state balance of 22× flux ÷ clearance.)")

# ============================================================================
# 3. Sensory-terminal read point (adventitia r_s = r₀ + 0.4–2 µm) and radial decay
# ============================================================================
# Read-point distance (2026-08-22 measurement-anchored, Luff 2000 [c167] PMID 10602088): source-b K⁺ is released from the
#   **VSMC surface**; the sensory terminal = **SP+ sensory C-fiber varicosity**, at **0.4–2 µm** from the VSMC
#   (most; SP- sympathetic varicosities closer, <1 µm; SP+ forms no junction with the SMC → µm-scale diffusion, not an nm cleft).
#   ——**not source-a's sympathetic varicosity↔sensory-terminal cleft of 20–150 nm [c49]** (nerve-nerve apposition).
#   Read points scan 0.4/1/2 µm, representative value 1 µm (median of the measured range).
print("\n[3] Terminal read-point (adventitia 0.4–2 µm, Luff 2000 [c167]) concentration")
for dgap_um in [0.4, 1.0, 2.0]:
    r_s = R_VESSEL + dgap_um * 1e-6
    u_t = cyl_u(r_s, J_AVG_MOL)
    u_p = cyl_u(r_s, J_STOC_flux)
    print(f"    {dgap_um:.1f} µm from SMC surface: time-average {u_t*1e3:.1f} µM | "
          f"STOC peak {u_p*1e3:.0f} µM")

# ============================================================================
# 3.5 Transient (STOC 100 ms) peak: cylindrical steady state is the lower bound, transient near-cluster is the upper bound
# ============================================================================
# The cylindrical steady-state solution assumes the flux persists ~λ²/D ≈ 100 s to reach steady state; but an STOC lasts only 100 ms, so K⁺ diffuses only
#   √(D·t) ≈ √(7.4e-10 × 0.1) ≈ 8.6 µm, not reaching λ=272 µm → the transient peak is above steady state.
# Transient point source (STOC-cluster sub-micron oasis 0.2–0.75 µm², radius 250–490 nm, ZhuGe 2002 PMID 12084772; near-cluster approximation): continuous point source running for 100 ms
#   u(r,t) = Q̇/(α·4π·D·r)·erfc(r/√(4Dt)), Q̇ = I_STOC/F = 2.59e-15 mol/s.
#   (no k't clearance term: 100 ms ≪ 1/k=100 s, k't=1e-3, correction <0.05%, see Nicholson Eq 10.2)
print("\n[3.5] Transient (STOC 100 ms) peak: cylindrical steady state (lower bound) vs transient point source (near-cluster upper bound)")
from scipy.special import erfc
Q_dot = I_STOC / F                                # mol/s point-source release rate
for dgap_um in [0.4, 1.0, 2.0]:
    rr = dgap_um * 1e-6                           # distance from STOC cluster, m
    arg = rr / np.sqrt(4.0 * D_eff * T_STOC)
    u_trans = Q_dot / (ALPHA * 4.0 * np.pi * D_eff * rr) * erfc(arg)
    u_cyl   = cyl_u(R_VESSEL + rr, J_STOC_flux)   # cylindrical steady-state comparison
    print(f"    {dgap_um:.1f} µm from STOC cluster: transient point source {u_trans*1e3:.0f} µM"
          f" (cylindrical steady-state lower bound {u_cyl*1e3:.0f} µM)")
print("    → True STOC peak lies between the two: ~0.10 mM (cylindrical uniform steady-state lower bound, 1 µm read point) – ")
print("      ~0.9 mM (1 µm) / ~2.2 mM (0.4 µm) near-cluster transient upper bound; STOC is a localized source (point-source treatment, ")
print("      S-independent) → near-cluster true value leans toward the transient point-source upper bound (0.4–1 µm read point ~0.9–2.2 mM); time-average ~10 µM.")

# ============================================================================
# 4. Density scaling (along-meridian 1.5× / acupoint 3.0×) + integrator (τ_m = 110 ms)
# ============================================================================
# Density (2026-08-23 correction, dropping the old 4.8×):
#   The old 4.8× = acupoint PCV appearance rate 87.5%/non-acupoint 18.1% (Álvarez-Prats 2019, PMID 31015854),
#   but appearance rate ≠ density (under Poisson, a 4.8 appearance-rate ratio ≈ 10× density ratio; under clustering <4.8, direction unknown),
#   and it only measured acupoint vs a 4cm control point, not the between-acupoint meridian segments. So it is split into two coefficients:
#   · along-meridian 1.5× = along-line average (acupoint 3.0 × acupoint fraction along the meridian ~0.15–0.33 + between-acupoint segments ~1.0)
#   · acupoint 3.0× = acupoint-local microvessel density (conservative placeholder; Liu 2014/Hsiu 2007 only "significantly higher" with no factor;
#     no difference among CV3/4/5/7 on the Ren meridian, PMC8099720 → no extra amplification in the meridian segment)
#   Steady-state tissue-average [K⁺] ∝ source density (source rate per volume) × 1/k → density × → [K⁺] ×.
#   PSC propagates along the meridian → main value uses along-meridian 1.5×; acupoint 3.0× is the acupoint-local upper bound (output for comparison only).
# Recruitment: recruitment is a varicosity firing event (source a) and does not directly multiply source b (the vessel); the "recruitment" analog for source b is
#   cAMP/cGMP ↑ spark frequency, but the quiet resting state of acupuncture/qigong has no strong NO/β drive, so conservatively not multiplied.
DENSITY     = 1.5     # along-meridian coefficient (main value, PSC propagates along the meridian)
DENS_POINT  = 3.0     # acupoint coefficient (acupoint-local upper bound, output for comparison only)

def integrator(dK_mM):
    """Membrane integrator: a 100 ms rectangular K⁺ pulse of amplitude dK → peak depolarization ΔV_peak.
    ΔV_peak = Nernst(dK) × (1 − e^{−T_STOC/τ_m}), T_STOC=100 ms, τ_m=110 ms."""
    dV_raw = NERNST * np.log((K_REST + dK_mM) / K_REST)
    return dV_raw * (1.0 - np.exp(-T_STOC / TAU_M))

print("\n[4] Density (along-meridian 1.5× / acupoint 3.0×) + terminal integrator (STOC 100 ms vs τ_m 110 ms → 60% integration)")
r_s = R_VESSEL + 1.0e-6           # take 1 µm representative read point (median of Luff 2000 [c167] measured 0.4–2 µm)
u_stoc_r3  = cyl_u(r_s, J_STOC_flux)
u_tonic_r3 = cyl_u(r_s, J_AVG_MOL)
u_stoc_d   = u_stoc_r3 * DENSITY
u_tonic_d  = u_tonic_r3 * DENSITY
u_stoc_p   = u_stoc_r3 * DENS_POINT
dV_stoc_0  = integrator(u_stoc_r3)      # u is already in mM (mol/m³ ≡ mM)
dV_stoc_d  = integrator(u_stoc_d)
dV_stoc_p  = integrator(u_stoc_p)
print(f"    Ordinary skin: STOC peak {u_stoc_r3*1e3:.0f} µM → post-integration {dV_stoc_0:.2f} mV")
print(f"              time-average {u_tonic_r3*1e3:.1f} µM")
print(f"    Along-meridian 1.5×: STOC peak {u_stoc_d*1e3:.0f} µM ({u_stoc_d:.2f} mM)"
      f"→ post-integration {dV_stoc_d:.2f} mV")
print(f"              time-average {u_tonic_d*1e3:.0f} µM")
print(f"    Acupoint 3.0×: STOC peak {u_stoc_p*1e3:.0f} µM ({u_stoc_p:.2f} mM)"
      f"→ post-integration {dV_stoc_p:.2f} mV (acupoint-local upper bound, not the along-meridian propagation main value)")
print(f"    Reach threshold: Δ[K⁺]o needs {DK_FULL:.1f} mM ({DV_TH:.1f} mV)")

# ============================================================================
# 5. Three-tier ladder closed loop (quantitative boundaries of quiescence / dam-break / resonance)
# ============================================================================
print("\n[5] Three-tier ladder closed loop (source-b perspective)")
print("    ┌ Quiescence (healthy default): source-b flux + density, STOC peak below threshold, still sub-full-threshold")
print(f"    │   Ordinary: STOC peak {u_stoc_r3*1e3:.0f} µM (integrated {dV_stoc_0:.2f} mV)")
print(f"    │   Along-meridian: STOC peak {u_stoc_d*1e3:.0f} µM ≈ {u_stoc_d:.2f} mM"
      f" (integrated {dV_stoc_d:.1f} mV)")
print(f"    │   → along-meridian STOC peak {u_stoc_d:.2f} mM is below the full threshold {DK_FULL:.1f} mM")
print("    ├ Dam-break (pathological breakdown/high-resistance reservoir): clearance shut off, k↓ ~1000× → λ↑ → steady state ↑ ~1000×")
print("    │   → time-average also crosses threshold (sustained plateau, not a 10% duty-cycle pulse)")
print("    └ Resonance (training/recruitment reservoir): spark frequency↑ (cAMP/cGMP) + whole-cell synchronization → duty cycle 10%→~100%")
print("        → STOC peak becomes a sustained plateau")
duty_base = F_STOC * T_STOC
print(f"    Time-average = peak × duty cycle ({duty_base:.2f}): ordinary {u_stoc_r3*1e3*duty_base:.0f} µM, "
      f"along-meridian {u_stoc_d*1e3*duty_base:.0f} µM; when duty cycle→1 it equals the peak itself")

# ============================================================================
# 6. Parameter sensitivity (honest bounds)
# ============================================================================
print("\n[6] Sensitivity (STOC amplitude ±, α 0.2/0.3, frequency 0.5/1/2 Hz)")
print(f"    {'STOC amp(pA)':>9} {'α':>4} {'along-meridian STOC peak(µM)':>14} {'integrated(mV)':>8} {'time-average(µM)':>12}")
for I_s in [125.0, 250.0, 500.0]:
    for alpha in [0.2, 0.3]:
        for f_s in [0.5, 1.0, 2.0]:
            J_st = (I_s * 1e-12) / (F * S_CELL)
            u_p = cyl_u(r_s, J_st, alpha) * DENSITY
            u_avg = u_p * (f_s * T_STOC)
            dV = integrator(u_p)
            print(f"    {I_s:9.0f} {alpha:4.1f} {u_p*1e3:14.0f} {dV:8.2f} {u_avg*1e3:12.0f}")
print("    (frequency only changes the time-average/duty cycle, not the STOC peak; α=0.2 is 1.5× higher than 0.3)")

# ============================================================================
# 7. Side by side with source a (varicosity): which is the primary source
# ============================================================================
print("\n[7] Source a (varicosity) vs source b (vascular BKCa): which is the primary source")
# Source a: varicosity AP K⁺ efflux Q_B (Keynes & Ritchie 1965 ~1 pmol/cm²/AP × varicosity area).
Q_B = 5.4e-20           # mol/AP single-varicosity K⁺ efflux (1.0 pmol/cm² × 5.4 µm²)
F_B = 1.0               # Hz varicosity baseline firing
J_a_avg = Q_B * F_B / (5.4e-12)   # mol/m²/s (single-varicosity surface area 5.4 µm²)
print(f"    Source-a varicosity time-average flux = {J_a_avg*1e8:.2f} pmol/cm²/s")
print(f"    Source-b vascular time-average flux = {J_AVG_MOL*1e8:.1f} pmol/cm²/s")
print(f"    → source b / source a = {J_AVG_MOL/J_a_avg:.0f}× (per membrane area); source b is the primary source")
print(f"    ⚠️ Geometry: source a sits in a 20–150 nm cleft (near field ∝ 1/r), source b in the 0.4–2 µm adventitia"
      f" (Luff 2000 [c167] measured) ——")
print(f"       so source b has large flux but is far, source a has small flux but is near; their contributions to the terminal concentration must be computed separately.")

# ============================================================================
# 8. Overall conclusion (citing the actually computed values)
# ============================================================================
gap_stoc = DK_FULL / u_stoc_d          # factor by which the along-meridian STOC peak is below the full firing threshold
gap_cons = DK_FULL / u_tonic_d         # factor by which the along-meridian time-average is below the full firing threshold
print("\n" + "=" * 86)
print("Overall conclusion")
print("=" * 86)
print(f"""
  ① Conceptual correction (overturning the old 8.6 µM):
     The old "source b = 8.6 µM" came from C·ΔV/F = 0.26 pmol/cm² (the net charge of one 25 mV hyperpolarization).
     That is the wrong treatment: K⁺ efflux is a **sustained circulating current** (each STOC counters the inward leak current).
     A single STOC moves {J_STOC_MOL*1e8:.0f} pmol/cm² ≈ {RATIO_STOC:.0f}× the capacitance lower bound;
     the time-average flux {J_AVG_MOL*1e8:.0f} pmol/cm²/s ≈ {RATIO_AVG:.0f}× the capacitance lower bound per second.

  ② But "flux ×22" ≠ "concentration ×22" — concentration is set by the "flux ÷ clearance" steady-state balance:
     · time-average steady state (cylindrical source, α=0.3): u(r₀) ≈ {u_tonic_r0*1e3:.0f} µM
     · STOC-period quasi-steady state (10× flux): u(r₀) ≈ {u_stoc_r0*1e3:.0f} µM
     · i.e. source-b concentration ~{u_tonic_r0*1e3/8.6:.0f}× (time-average) / ~{u_stoc_r0*1e3/8.6:.0f}×
       (STOC peak) of the old 8.6 µM —— flux is 22×, but concentration only amplifies ~{u_stoc_r0*1e3/8.6:.0f}×,
       because clearance pulls the steady state back to the µM–sub-mM range.

  ③ Three-tier ladder closed loop (source-b perspective):
     · Ordinary skin: STOC peak {u_stoc_r3*1e3:.0f} µM, time-average {u_tonic_r3*1e3:.0f} µM
       → terminal integration {dV_stoc_0:.2f} mV, ~{DK_FULL/(u_stoc_r3):.0f}× below the full threshold, quiescent.
     · Along-meridian density 1.5× (acupoint 3.0×): STOC peak {u_stoc_d*1e3:.0f} µM ≈ {u_stoc_d:.2f} mM (below the full threshold {DK_FULL:.1f} mM)
       → terminal integration {dV_stoc_d:.1f} mV; time-average {u_tonic_d*1e3:.0f} µM (still {gap_cons:.0f}× below the full threshold).
     · Conclusion: density compresses the "ordinary→full threshold" gap from the old ~125× to ~{gap_stoc:.0f}×.
       The STOC peak is now **below the full threshold ({DK_FULL:.1f} mM)**; the time-average is still ~1 order of magnitude low.

  ④ The correct role of source b (repositioned):
     · Source b is the **primary source** (unit membrane-area flux = {J_AVG_MOL/J_a_avg:.0f}× that of source a),
       but geometrically far from the terminal (adventitia 0.4–2 µm) vs source a near (cleft <150 nm).
     · Source b's role = **raising the quiescent baseline from µM into the sub-threshold range** (along-meridian STOC peak {u_stoc_d:.2f} mM, uniform treatment),
       so that a "reservoir event" (dam-break/resonance) needs only ~{gap_stoc:.0f}× more (instead of 125×) to reach a sustained plateau.
     · Crossing threshold (sustained mM plateau) still relies on reservoir events: dam-break (clearance shut off, k↓~1000× → time-average also crosses threshold)
       or resonance (spark frequency↑/synchronization → duty cycle 10%→100% → STOC peak becomes a sustained plateau).
     · This "source b = sub-threshold charging, reservoir = closing the gate" division of labor is fully self-consistent with the three-tier ladder quiescence/dam-break/resonance.

  ⑤ Honest bounds:
     · The 250 pA STOC amplitude (Benham & Bolton −40 mV) has reduced driving force at the relaxed V_m≈−55 mV,
       so the actual flux may be ~2–5× lower (sensitivity table section [6]); α=0.2 gives 1.5× higher concentration.
     · Recruitment (~2×) is a varicosity firing event of source a, not multiplied directly onto source b; if "spark frequency↑" is treated as recruitment,
       it affects the time-average/duty cycle (toward resonance), not the STOC peak.
     · Robustness of the conclusion: whatever the STOC amplitude/α, source b + density lands on the "sub-threshold mM edge",
       not crossing the sustained threshold; the qualitative conclusion (source b primary source + sub-threshold charging + reservoir closing the gate) does not depend on specific parameters.
""")
