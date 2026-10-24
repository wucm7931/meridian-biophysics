# SI Table S3 — Hyperpolarized-type molecular subtypes: mechanistic details and primary literature

> **Related sections**: main text Section III C (hyperpolarized-type overview), Section III F (hyperpolarized-type malignant pseudo-A+ pathway)
> **Classification**: H1=Ca²⁺-driven, H2=neural-crest/neuroendocrine-origin, H3=K⁺-channel-overexpression, H4=HIF-driven
> **Symbols**: R = \|RMP\|/100 (defined in SI Table S1); R_C/R_N = \|RMP\|_C/\|RMP\|_N (the factor 100 cancels) = cancer/normal RMP ratio. The "RMP range (mV)" column of the second table is the measured/reported RMP range (signed value) of the cells in that row, complementary to the \|RMP\|_N/\|RMP\|_C (mV) absolute-value convention of the first table.

---

## Hyperpolarized-type confirmed case (1 pair, ratio > 1.0)

| Subtype | Cancer | Cell line | \|RMP\|_N (mV) | \|RMP\|_C (mV) | R_C/R_N | Core molecular mechanism | PMID | Source |
|:---:|------|------|:---:|:---:|:---:|------|------|------|
| H1 | colon cancer | HT-29 | 45 | 64 | 1.42 | ER Ca²⁺ depletion → TRPC1↑/STIM2↓ (reciprocal transfer) + ORAI1↑/STIM1↑ → SOCE↑ (ICRAC+ISOC mixed current) → Ca²⁺ influx↑ → hyperpolarization (Ca²⁺-activated K⁺ conductance; KCa3.1 is the downstream inferred channel); hyperpolarization enlarges the Ca²⁺ driving force, forming a positive-feedback loop | 25143380 | Sobradillo 2014 J Biol Chem |

> ⚠️ **The excluded "hyperpolarized" candidate — medullary thyroid carcinoma (MTC)**: MTC (normal 36.6 → cancer 50.1 mV, ratio 1.37) was originally listed as an H2 confirmed case, and is now excluded from the clinical hyperpolarized cases per SI Table S1 Note 8 — its "hyperpolarization" is a cross-lineage + cross-species artifact: the normal value 36.6 mV comes from the hamster thyroid follicular epithelium (Jamakosmanovic & Loewenstein 1968), and the cancer value 50.1 mV from the human thyroid parafollicular C cell (neuroendocrine origin, Sand 1986); the two are not the same cell lineage, and the C cell's polarized electrical phenotype is a developmental retention rather than cancerous hyperpolarization. The H2 subtype (neural-crest/neuroendocrine-origin) is now anchored only on the SCLC neuroendocrine cell lines (see the "dynamic/experimental cases" table below).

## Hyperpolarized-type dynamic/experimental cases (switchable ratio)

| Subtype | Cancer | Cell line | RMP range (mV) | Switching mechanism | PMID | Source |
|:---:|------|------|:---:|------|------|------|
| H1→H4 | melanoma | IGR1 | −22.3 ~ −77.0 | Ca²⁺ present → KCa open → −77 mV (ratio=1.10); Ca²⁺ removed → KCa closed → −22.3 mV (ratio=0.32). A single cell line spans a dynamic range of 55 mV (note: the closed state −22.3 mV is from Gavrilova-Ruch 2002; Nilius 1990 reports a resting value of about −30 mV, which is the value used in SI Table S1 row 25) | 12172639; 2177183 | Gavrilova-Ruch 2002; Nilius 1990 |
| H4 | melanoma | IGR1 (HIF) | −77 (HIF→KCa↑) | HIF-1α → SK2/IK transcriptional upregulation → hyperpolarization; VHL overexpression → HIF-1α degradation → KCa↓ → depolarization reversal | 17065157 | Tajima 2006 |
| H2 | SCLC (NE) | NE cell lines | −72 (NE phenotype retained, no normal reference) | neuron-like electrical excitability: Nav+Kv+KCa → action potentials; electrical activity causally drives malignancy (ChR2 → depolarization → tumor doubling) | 39939778 | Peinado 2025 Nature |
| H2 | SCLC (brain metastasis) | 16T | −72 (GABA reversal −27.3) | brain microenvironment → GABA depolarization (−27.3 mV reversal potential) → altered Cl⁻ homeostasis; neuron–cancer-cell synapse | 40931074 | Savchuk 2025 Nature |
| H3 | TNBC | MDA-MB-231 (K⁺ OE) | baseline −30 → hyperpolarized (ratio 0.67→>1.0) | K⁺-channel overexpression → hyperpolarization → CDH11↑ → MAPK → migration/metastasis↑. Demonstrates that the depolarized/hyperpolarized boundary can be crossed | 34954628 | Payne 2021 |

> ⚠️ **Table notes (H2/H3 reference convention)**: ① **The two H2 rows (SCLC) set no normal reference** — SCLC originates from pulmonary neuroendocrine cells (PNEC), and its −72 mV is a **developmental polarized-phenotype retention** of the neuroendocrine lineage (PNEC neuron-like electrical excitability), not cancerous hyperpolarization "relative to normal alveoli"; the "normal alveolar −50" reference has no single traceable source and contradicts the measurements (SI Table S1 Note 9: rat AT2 measured −4.4 mV, naturally depolarized). This is the same logic as the cross-lineage exclusion of MTC in Note 8. ② **The denominator of the H3 ratio is uniformly the normal mammary 45 mV** (S1 Note 1): the MDA-MB-231 baseline adopts Quicke −30 (ratio 0.67, S1 Note 7), and Payne −19.5 is the alternative measurement (ratio 0.43). The old "0.72 = −19.5 / MCF10A −27.1" took Payne's internal MCF10A as the denominator, inconsistent with the 45 mV convention of this table, and has been abandoned.

---

> **References** (SI-S3)
> 1. Sobradillo D et al. J Biol Chem 2014;289(42):28765–28782. PMID: 25143380
> 2. Sand O et al. Acta Physiol Scand 1986;126(2):173–179. PMID: 3705981
> 3. Gavrilova-Ruch O et al. J Membr Biol 2002;188:137. PMID: 12172639
> 4. Nilius B et al. Pflügers Arch 1990;417:269–277. PMID: 2177183
> 5. Tajima N et al. J Biol Chem 2006;281(52):39925–39934. PMID: 17065157
> 6. Peinado P et al. Nature 2025;638(8052):1102–1110. PMID: 39939778
> 7. Savchuk S, Gentry KM, et al. Neuronal activity-dependent mechanisms of small cell lung cancer pathogenesis. Nature 2025;646(8087):1232–1242. PMID: 40931074
> 8. Payne SL et al. EBioMedicine 2021;75:103767. PMID: 34954628
