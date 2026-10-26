# SI Data 4: Same-lineage RMP/R gradients: primary-literature data excerpts

> **Related section**: main text Section III C (R and cell lifespan: the monotonic gradient within a lineage)
> **Principle**: across types, \|RMP\| has no statistical relation to lifespan; within a single lineage, R increases monotonically along the differentiation/maturation axis
> **Data status**: ⚠️ same-lineage data points are sparse and measurement methods are not standardized. The following are the currently available data; the conclusion belongs to the mechanistic-inference layer with experimental anchoring.
> **Symbols**: R = \|RMP\|/100 (defined in Data 1).

---

## Lineage 1: mesenchymal stem cell → osteogenesis (same lineage)

| State | R | RMP (mV) | Method | Channel mechanism | PMID | Source |
|------|:---:|:---:|------|------|------|------|
| Undifferentiated MSC (OS baseline, week 0) | 0.37 | −37.0±9.4 | intracellular recording | unidentified (does not specify the channel) | 19011685 | Sundelacruz 2008 PLoS ONE |
| Osteogenic differentiation week 4 | 0.93 | −93.0±5.5 | DiSBAC2(3) dye extrapolation | unidentified | 19011685 | ibid. |
| Intervention: depolarization (ouabain / high extracellular K⁺) | — | — | — | depolarization → differentiation inhibition | 19011685 | **Blocking hyperpolarization = blocking differentiation** (causal evidence) |

> ⚠️ Osteogenesis is a same-lineage (intra-mesodermal) maturation; the −37/−93 mV are from Sundelacruz 2008 (week 0 is intracellular recording; weeks 1–4 are DiSBAC2(3) dye extrapolated at "1% fluorescence ≈ 1 mV", so the absolute-value precision is limited, and the AD branch already crosses E_K≈−90, so only the OS branch is used).

## Lineage 2: BMSC → hair-cell-like (cross-germ-layer transdifferentiation)

| State | R | RMP (mV) | Method | Channel mechanism | PMID | Source |
|------|:---:|:---:|------|------|------|------|
| BMSC → hair-cell-like transdifferentiation D17 | 0.12 | −11.9±0.9 | whole-cell patch clamp | — | 41016575 | Peyvandi 2025 Exp Cell Res |
| BMSC → hair-cell-like transdifferentiation D21–26 | 0.59 | −59.0±1.1 | whole-cell patch clamp | — | 41016575 | ibid. |

> ⚠️ The hair-cell-like is a **cross-germ-layer transdifferentiation** (mesodermal MSC → ectodermal hair cell), not a "proliferation→differentiation" developmental axis within a single normal lineage; the −11.9/−59.0 mV are from Peyvandi 2025 whole-cell patch clamp (plasma-membrane RMP, not mitochondrial). Its value lies in that, even across germ-layer transdifferentiation, the differentiation direction remains consistent with the same-lineage case (differentiation → hyperpolarization → R↑), serving only as a generalization corroboration of "the differentiation direction is consistent", not a test of the "same-lineage R↑" hypothesis.

## Lineage 3: small-intestinal epithelium: crypt → villus-tip differentiation gradient

| State | R | RMP (mV) | Method | Channel mechanism | PMID | Source |
|------|:---:|:---:|------|------|------|------|
| Crypt stem cell | 0.15–0.40 | −15~−40 | microelectrode | K⁺ channels low / undifferentiated | 6149040; 2782407 | multi-source summary (crypt–villus heterogeneity) |
| Villus-tip absorptive cell | 0.70–0.72 | −70~−72 | microelectrode | Kir/K2P upregulation | 6149040; 2782407 | multi-source summary |

⚠️ The continuous RMP-gradient data along the intestinal crypt–villus axis come from multiple animal models (rat, rabbit, guinea pig); human in-situ data are extremely scarce.

## Lineage 4: epidermis: basal → stratum corneum differentiation

| State | R | RMP (mV) | Method | Channel mechanism | PMID | Source |
|------|:---:|:---:|------|------|------|------|
| Basal layer (proliferating) | 0.24–0.27 | −24~−27 | — | low K⁺ conductance | 1690740; 17359336 | field reference value (Mauro 1990; Gönczi 2007) |
| Spinous/granular layer (differentiating) | 0.42 | −42 | — | Kir/K2P progressive upregulation | 11149902 | field reference value (Koegel 2001) |

⚠️ The methodological challenge of epidermal RMP measurement is enormous (keratinized barrier, mixed cell layers); the existing values are mainly animal-model data and in-vitro keratinocyte cultures.

## Lineage 5: glial cells: glioma → normal glia (⚠️ cancer–normal gradient, not a developmental lineage)

> ⚠️ **Nature of the data**: this "lineage" is actually the **glioma grades I–IV vs. normal glia** measured by Picker 1981 within the same study — a cancer–normal gradient, not a "proliferation→differentiation" developmental axis within a single normal glial lineage. It overlaps completely with the glioma-grade gradient of main text Section III C, and is retained here only as cross-lineage directional corroboration, not as independent developmental-lineage evidence.

| State | R | RMP (mV) | Method | Channel mechanism | PMID | Source |
|------|:---:|:---:|------|------|------|------|
| Glioblastoma (grade IV) | ~0.40 | ~−40 | microelectrode | low Kir4.1 | 6267226 | Picker 1981 J Neurosurg |
| Normal glia (reference) | ≈0.69 | ≈−69 | microelectrode | high Kir4.1 expression | 6267226 | ibid. |

## Lineage 6 (reference): neurons: postnatal development (infant → young adult)

| State | R | RMP (mV) | Method | Channel mechanism | PMID | Source |
|------|:---:|:---:|------|------|------|------|
| Neonatal (infant <1 y) | 0.61 | −60.6 | patch clamp | low K⁺ leak conductance | 40152903 | Barzó 2025 eLife |
| Mature (young adult 20–39 y) | 0.69 | −68.7 | patch clamp | K⁺ leak conductance↑ (KCNK) | 40152903 | ibid. |

---

## Methodological notes

1. **Cross-lineage incomparability**: the absolute values, measurement techniques, and species differ across lineages, so the absolute levels of R cannot be compared directly.
2. **Directional consistency**: in the five lineages with available data (plus one cross-lineage transdifferentiation corroboration), R increases monotonically with differentiation/maturation (lineage 2 is cross-germ-layer transdifferentiation — directionally consistent but not a developmental axis within a single normal lineage).
3. **Causality**: only the MSC lineage has direct intervention experiments (depolarization [ouabain / high K⁺] → differentiation inhibition, hyperpolarization [pinacidil / diazoxide] → differentiation promotion; Sundelacruz et al., 2008); the other lineages are correlational.
4. **Data gaps**: the hematopoietic, hepatic, and pancreatic lineages lack quantitative RMP differentiation-gradient data.

---

> **References** (SI Data 4)
> 1. Sundelacruz S, Levin M, Kaplan DL. PLoS ONE 2008;3(11):e3737. PMID: 19011685 DOI: 10.1371/journal.pone.0003737 (= main reference [119])
> 2. Picker S et al. J Neurosurg 1981;55:347–363. PMID: 6267226 DOI: 10.3171/jns.1981.55.3.0347 (= main reference [116])
> 3. Barzó P, Szöts I, Tóth M, et al. eLife 2025;13:RP100390. PMID: 40152903 DOI: 10.7554/eLife.100390 (= main reference [110])
> 4. Mauro TM, et al. 1990 (membrane potential and differentiation of human keratinocytes cultured in low calcium). PMID: 1690740 DOI: 10.1002/jcp.1041430103
> 5. Koegel H, Alzheimer C. Expression and biological significance of Ca²⁺-activated ion channels in human keratinocytes. FASEB J 2001. PMID: 11149902 DOI: 10.1096/fj.00-0055com
> 6. Gönczi M, et al. Exp Dermatol 2007;16:302–310. PMID: 17359336 DOI: 10.1111/j.1600-0625.2006.00533.x
> 7. Positional dependence of enterocyte membrane potential in hamster and rabbit enterocytes. Comp Biochem Physiol 1984. PMID: 6149040 DOI: 10.1016/0300-9629(84)90613-3
> 8. A microelectrode study of responses to secretagogues by epithelial cells on villus and crypt of rat small intestine. Am J Physiol 1989. PMID: 2782407 DOI: 10.1152/ajpgi.1989.257.3.G334
> 9. Peyvandi AA, Davoudi S, Bazgir N, et al. A time course analysis of the electrophysiological and gene expression properties during differentiation of hair cell-like cells in culture. Exp Cell Res 2025;452(2):114775. PMID: 41016575 DOI: 10.1016/j.yexcr.2025.114775 (= main reference [252])
