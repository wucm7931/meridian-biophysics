# -*- coding: utf-8 -*-
"""
model_v2 · shared physical constants and tissue parameters (single source of truth)
=====================================================
Design iron rules (see 00_overview.md §3):
  · Each numeric literal is annotated on the same line with "source (PMID/local file) + original measured value".
  · Derived values (class B) are computed on the spot from the original values; do not hardcode the result.
  · Unit rule: 1 mol/m³ ≡ 1 mM; 1 pmol/cm² = 1e-8 mol/m²; 1 µm² = 1e-12 m².

Physical constants (F, R, T) are SI defined quantities, not literature measurements; the rest are literature measurements (class A)
or quantities derived from class A (class B). Imported by each submodule script.
"""

# ----------------------------------------------------------------------
# Physical constants (SI defined quantities, no literature source)
# ----------------------------------------------------------------------
F      = 96485.0          # Faraday constant [C/mol]
R_gas  = 8.314            # Gas constant [J/(mol·K)]
T      = 310.0            # Physiological body temperature 37°C [K]
N_A    = 6.02214076e23    # Avogadro constant [/mol] (SI defined quantity, for converting ion counts)

# ----------------------------------------------------------------------
# Class B · Nernst coefficient (derived from R/T/F)
#   Formula: dV = (R·T/F)·ln([K]2/[K]1); ×1000 converts V→mV
# ----------------------------------------------------------------------
NERNST = R_gas * T / F * 1000.0     # = 26.71 mV (RT/F natural-log slope, 37 °C)

# ----------------------------------------------------------------------
# Class A · free-water KCl diffusion coefficient
#   Original measurement: D_KCl = 1.9×10⁻⁵ cm²/s (Syková & Nicholson 2008, Physiol Rev
#   88:1277–1340, PMID 18923183)
# ----------------------------------------------------------------------
D_free = 1.9e-9            # m²/s (= 1.9e-5 cm²/s)

# ----------------------------------------------------------------------
# Class A · ECS tortuosity
#   Original measurement: λ = 1.6 (brain ECS; Nicholson & Syková 1998, Trends Neurosci
#   21:207–215, PMID 9610885)
# ----------------------------------------------------------------------
lam_T  = 1.6

# ----------------------------------------------------------------------
# Class B · effective diffusion coefficient (derived D_eff = D_free / λ²)
# ----------------------------------------------------------------------
D_eff  = D_free / lam_T**2     # = 7.42e-10 m²/s = 742 µm²/s

# ----------------------------------------------------------------------
# Class C · ECS volume fraction (local porosity at the source b read point) — conservative lower bound 0.30
#   Physics: α is the diffusible-fluid volume fraction at a read point 0.4–2 µm outside the arteriole
#   (inside the perivascular interstitial cuff), not the bulk tissue porosity. u_b ∝ 1/α, so smaller α → larger u_b.
#   Three-layer geometric evidence (local full text, see 02 §6):
#     ① Negrini 2001, PMID 11704543 (rabbit lung): perivascular interstitial cuff around ~30 µm arterioles has mean thickness
#        25.1±7.1 µm (≈1.7× vessel radius) → read point lies inside the open cuff, not in brain-parenchyma packing.
#     ② Hahn 2025, PMID 41085930 (skin): interstitial volume = 50% of skin mass (muscle only 10%)
#        → skin interstitium is open; brain ECS α=0.2–0.24 is an inappropriate analogy.
#     ③ Benias 2018, PMID 29588511 (human, including dermis + vascular adventitia): the adventitia is not a dense collagen wall
#        but a fluid-filled interstitial space (collagen bundles lined with cells on one side, exposed to fluid on the other) → local α direction > 0.3.
#   Conclusion: the three papers give only structural geometry and none gives a dimensionless porosity → the true α > 0.30 but is not directly measured;
#   take 0.30 as the conservative lower bound (the upper value for brain ECS): α↑ → u_b↓, so the lower bound overestimates source b;
#   if u_b(α=0.30) is still < 6.8 mM the conclusion is robust. See each submodule for sensitivity.
# ----------------------------------------------------------------------
ALPHA  = 0.30

# ----------------------------------------------------------------------
# Class C · first-order clearance rate constant (placeholder: main value 0.01 /s, range [0.003, 0.04] /s)
#   Anchoring: in-vivo muscle recovery τ≈270–360 s (Nielsen 2003, PMID 14634198) → k≈0.003–0.004 /s;
#         muscle-pump net clearance tangent dR/d[K+] = Vmax·Km/([K+]+Km)² = 0.09×3.8/7.7² ≈ 0.006 /s
#         (Lo 2004, PMID 15277205); spinal cord normal 4.6e-3 /s (Syková 1994, PMID 8113325).
#   Upper-end reference: in-vitro SCG electrical-stimulation K⁺ decay τ≈25 s → k = 0.04 /s
#         (Förstl et al. 1982, Neuroscience 7:3221–3229, PMID 7162635).
#   The headline predictions (quiescence / independent charging / speed) are robust to k; see memo_k_clearance_muscle.md.
# ----------------------------------------------------------------------
K_CLEAR = 0.01               # /s (main value; range 0.003–0.04)

# ----------------------------------------------------------------------
# Class A · resting extracellular K⁺
#   Interstitial resting [K⁺]o = 3.5–5 mM (Bennet 2021 [c321] Clin Transl Med 11:e461, human interstitial fluid
#   ISF K⁺ range, taking the median 4.0 mM). Corroborating: MacLean 2000 PMID 10712273.
# ----------------------------------------------------------------------
K_REST  = 4.0                # mM

# ----------------------------------------------------------------------
# Class A · specific membrane capacitance
#   Standard physiology: C_m = 1 µF/cm² (textbook value, no single PMID)
# ----------------------------------------------------------------------
C_M     = 1e-2               # F/m² (= 1 µF/cm²)
