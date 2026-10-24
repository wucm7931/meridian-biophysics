# Derived number → script mapping

| Derived number | Main text/appendix location | Script |
|---|---|---|
| λ = √(D_eff/k) ≈ 272 µm (steady-state reaction–diffusion decay length) | Appendix 1.3.1, Appendix 4.2 | Constant source `common.py` (`D_eff=7.42e-10 m²/s`, `K_CLEAR=0.01`); computed in `model_source_a.py`, `model_source_b.py`, `model_arteriole_k.py` |
| Source-a varicosity S_v / Q_b / linear density; source-a bulk steady state ≈ 5–42 nM | Table 2.3, Appendix 1.3.2 | `model_sympathetic.py` (S_v, Q_b), `model_source_a.py` (cylindrical line-source K₀ solution u_bulk) |
| Source-b STOC (250 pA × 100 ms), main-source determination | Table 2.3, Appendix 1.3 | `model_arteriole_k.py`, `model_source_b.py` |
| Near-field transient peak ~1–9 µM (u_a_mid ≈ 3.4 µM), bulk steady state ~6 nM | Table 2.3, Appendix 1.3.2 | `calc_nearfield.py` |
| Bifurcation Λ\*≈15.89, S\*(Λ) saddle-node bistability, D(p) shape | Appendix 1.4.6 | `calc_bifurcation.py` (γ=0.4/n=2/P_eff=16) |
| Five-domain relay table (S converges from seed 0.8 toward S\*) | Appendix 1.4.5 | `calc_relay_table.py` (γ=0.4/n=2/P_eff=16/a=0.40) |
| Form robustness (18 combinations d×p×D) | Appendix 1.4.6.2 | `test_form_robustness.py` (7-equation framing γ=0.5; conclusion framing-independent) |
| Parameter robustness (12-parameter sweep) | Appendix 1.4.7 | `test_param_robustness.py` (7-equation framing γ=0.5; conclusion framing-independent) |
| D(p) power-law form, 8 shapes (pow1/2/3/half/pure/log_thresh…) | Appendix 1.4.6.2 | `test_power_form.py` (γ=0.5; form conclusion framing-independent, n<1 already re-checked at γ=0.4) |
| Health-layer four-segment profile (h0=10/γ=0.4/n=2/P_max=16) | Appendix 1.4.12 | `calc_onset_data.py` |
| Dam-break / uniform-room contrast | Appendix 1.3, §2.3.5 | `model_uniform_room.py` |
| Configurational-entropy curve S_config/k_B, \|E_K\| curve, three example states (f_i=0.8→\|RMP\| 0, 0.947→~40 mV, 0.981→~69 mV) | Figure 2.4-1, Appendix 2.1.8 | `fig_entropy_rmp.py` (computed figure, α=0.5/0.2 contrast) |
| 33 cancer/normal R pairs, sign test, five p-values (3.0×10⁻⁸ / 2.4×10⁻⁴ / 4.0×10⁻⁵ / 3.0×10⁻⁶ / 7.4×10⁻³) + counts | 04 §2.4.2, Appendix 2.3, Appendix 2.4 | `stat_rmp_cancer.py` (reads `data/paired_rmp_data.csv` + `data/srivastava_2021_extracted.csv`, two-sided exact binomial) |
| Five n-value table (Λ\*, S\*(Λ=20), h_low/h_crit three columns) | Appendix 1.4.7 | Λ\* column `calc_bifurcation.py`; S\*(Λ=20) column `calc_n_sweep.py`; h_low/h_crit column `calc_n_window.py` (P_max=16 final framing) |
| Critical-slowing prefactor C=π/√(H_Λ\|H_SS\|/2)≈10.96 | Appendix 1.4.6.3 | `calc_critical_slowdown.py` |
| Short-range decay segment counts (gap→segments: 8→1 / 1→9 / 0.2→22 / 0.05→47 / 0.01→108) | Appendix 1.4.12.2 | `calc_decay_segments.py` |
