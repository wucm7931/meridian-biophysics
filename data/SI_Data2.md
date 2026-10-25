# SI Data 2: Extraction and organization of the paired data from the Srivastava 2021 meta-analysis

> **Related sections**: main text Section III C (Srivastava regression), Methods Section VII B (regression extraction)
> **Data source**: Srivastava P, Kane A, Harrison C, Levin M. *Bioelectricity* 2021;3(1):42–67. PMID: 34476377 (Supplementary Table S1)
> **Extraction method**: extracted row by row directly from the original Supplementary Table S1, with cross-checking against the primary literature (see the QT-35 correction in Note 3)

---

## Paired data (for regression validation)

The mixed-effects regression reported by Srivastava 2021 is $V_{mem}^{cancer} - V_{mem}^{somatic} = -16.9 - 0.85 \cdot V_{mem}^{somatic}$ ($p < 0.0001$, 41 cancer types, 70 normal cell types). This table contains the core rows with normal pairing, **expanded row by row** (the 5 glioma grades and the 11 rat-thyroid-cancer lines are listed as independent rows), for validating the consistency of the regression with this paper's R metric (Eqs. (5)(6)).

| Pair ID | Species | Tissue | Normal | \|RMP\|_N (mV) | Cancer | \|RMP\|_C (mV) | Ratio | Source |
|:---:|------|------|------|:---:|------|:---:|:---:|------|
| CX1 | human | cervix | endocervix | 33 | cervical cancer | 22 | 0.67 | Lash 1955 |
| GA1 | human | stomach | gastric-fundus mucosa | 56 | gastric cancer | 24 | 0.43 | Tokuoka 1957 |
| TH1 | human | thyroid | normal thyroid | 36.6 | medullary thyroid carcinoma | 50.1 | **1.37 (hyperpolarized)** | Jamakosmanovic 1968 / Sand 1986 |
| BR1 | human | breast | MCF-10A | 58.1 | MCF-7 | 42.1 | 0.72 | Marino 1994 |
| BR2 | human | breast | MCF-10A | 58.1 | MDA-435 | 51.8 | 0.89 | Marino 1994 |
| BR3 | human | breast | MCF-10A | 58.1 | intraductal carcinoma | 13.3 | 0.23 | Marino 1994 |
| GL1 | human | brain | normal glia | 69 | astrocytoma grade I | 70 | 1.01 (tied) | Picker 1981 |
| GL2 | human | brain | normal glia | 69 | astrocytoma grade II | 59 | 0.86 | Picker 1981 |
| GL3a | human | brain | normal glia | 69 | mixed glioma grade III (astrocytic) | 51 | 0.74 | Picker 1981 |
| GL3o | human | brain | normal glia | 69 | mixed glioma grade III (oligodendroglial) | 53 | 0.77 | Picker 1981 |
| GL4 | human | brain | normal glia | 69 | glioblastoma grade IV | 40 | 0.58 | Picker 1981 |
| SK1 | rodent | skin | mouse skin | 31.5 | squamous-cell carcinoma | 18.5 | 0.59 | Chowdhury 1973 |
| TH2 | rodent | thyroid | hamster thyroid | 38.5 | thyroid cancer | 22.6 | 0.59 | Jamakosmanovic 1968 |
| LI1 | rodent | liver | rat hepatocyte | 37.1 | liver cancer | 19.8 | 0.53 | Binggeli 1980 |
| FI1 | rodent | fibroblast | mouse corneal fibroblast | 42.5 | fibrosarcoma | 14.3 | 0.34 | Binggeli 1980 |
| FI2 | rodent | fibroblast | quail fibroblast | 40.1 | QT-35 fibrosarcoma | 20.7 | 0.52 | Binggeli 1985 (corrected, see Note 3) |
| TH3-1 | rodent | thyroid | rat thyroid | 47.2 | thyroid cancer 1-1C2 | 31.7 | 0.67 | Jamakosmanovic 1968 |
| TH3-2 | rodent | thyroid | rat thyroid | 47.2 | thyroid cancer 1-1D | 26.8 | 0.57 | Jamakosmanovic 1968 |
| TH3-3 | rodent | thyroid | rat thyroid | 47.2 | thyroid cancer 1-8 | 26.4 | 0.56 | Jamakosmanovic 1968 |
| TH3-4 | rodent | thyroid | rat thyroid | 47.2 | thyroid cancer 1-4 | 26.6 | 0.56 | Jamakosmanovic 1968 |
| TH3-5 | rodent | thyroid | rat thyroid | 47.2 | thyroid cancer 1-5A | 24.4 | 0.52 | Jamakosmanovic 1968 |
| TH3-6 | rodent | thyroid | rat thyroid | 47.2 | thyroid cancer 1-5E | 23.6 | 0.50 | Jamakosmanovic 1968 |
| TH3-7 | rodent | thyroid | rat thyroid | 47.2 | thyroid cancer 1-6 | 31.7 | 0.67 | Jamakosmanovic 1968 |
| TH3-8 | rodent | thyroid | rat thyroid | 47.2 | thyroid cancer 1-7 | 21.7 | 0.46 | Jamakosmanovic 1968 |
| TH3-9 | rodent | thyroid | rat thyroid | 47.2 | thyroid cancer 16-1 | 23.8 | 0.50 | Jamakosmanovic 1968 |
| TH3-10 | rodent | thyroid | rat thyroid | 47.2 | thyroid cancer 16-4 | 25.9 | 0.55 | Jamakosmanovic 1968 |
| TH3-11 | rodent | thyroid | rat thyroid | 47.2 | thyroid cancer 1-3 | 26.6 | 0.56 | Jamakosmanovic 1968 |
| MU1 | rodent | muscle | rat muscle fiber | 95.6 | rhabdomyosarcoma | 16.3 | 0.17 | Balitsky 1964 |
| 3T1 | rodent | 3T3 | Swiss 3T3 | 12.3 | Py3T3 transformant | 35.2 | **2.86 (hyperpolarized)** | Killion 1984 |

> **Note**: Srivastava 2021's original Supplementary Table S1 recorded the QT-35 cancer value as −40.1 mV (identical to the normal, ratio=1.00). After checking the original paper Binggeli & Weinstein 1985 (PMID 3965134), the correct QT-35 value under high-density conditions is −20.7 mV, and this table adopts the corrected value (see Data 1 Note 3 for details). The 2 "hyperpolarized" pairs in the table (TH1 medullary thyroid carcinoma, 3T1 Py3T3) are **not counted as clinical hyperpolarized cases**: TH1 is a cross-lineage + cross-species artifact (normal = hamster thyroid follicular epithelium vs. cancer = human parafollicular C cell), and 3T1 is an SV40 experimental transformation (the 3T3 baseline −12.3 mV is already an immortalized cell, not normal physiology) — the clinical hyperpolarized case is only HT-29 (see Data 1).

---

## Sign test (two conventions)

> All sign-test p values in this table (conventions A/B and sensitivity analysis) were recomputed with a two-sided exact binomial test (`scipy.stats.binomtest`).

### Convention A: row-by-row expansion (n = 29 pairs)

The 5 glioma grades and the 11 rat-thyroid-cancer lines are each treated as independent pairs. Of the 29 pairs, **26 have cancer R < normal R, 2 have cancer R > normal R (TH1, 3T1), and 1 is tied (GL1 grade I, \|ΔR\| = 0.01)**. Non-tied pairs n = 28, two-sided exact binomial test: depolarized 26/28 = 92.9%, **p = 3.0 × 10⁻⁶**.

### Convention B: same-normal-value pairs merged by mean (n = 15 pairs)

The 5 grades of GL_all are merged (mean 54.6 mV), and the 11 lines of TH3 are merged (mean 26.3 mV). Of the 15 pairs, **13 depolarized, 2 hyperpolarized (TH1, 3T1), 0 tied**. Non-tied pairs n = 15, two-sided exact binomial test: depolarized 13/15 = 86.7%, **p = 7.4 × 10⁻³**.

### Sensitivity analysis: excluding the 2 disputed hyperpolarized pairs (TH1, 3T1)

| Convention | Valid pairs | Depolarized | Tied | Non-tied | p |
|:---:|:---:|:---:|:---:|:---:|:---:|
| A expanded | 27 | 26 | 1 | 26 | 3.0 × 10⁻⁸ |
| B merged | 13 | 13 | 0 | 13 | 2.4 × 10⁻⁴ |

### Notes on the choice of convention

- **Convention A (expanded)**: retains all grade/line information (the glioma grade I→IV gradient is the core of the "graded depolarization" evidence in the main text), with larger n and smaller p; but the 5 GL grades share the same normal value −69, and the 11 TH3 lines share the same normal value −47.2, constituting partial pseudoreplication of "a shared normal reference within the same study", not strictly 29 independent measurements.
- **Convention B (merged)**: one pair per study per tissue, eliminating the pseudoreplication, more conservative; but it erases the glioma-grade gradient, with smaller n and larger p.
- The main-text statistics of this paper adopt the 33 pairs of Data 1 (same-lineage, more rigorously traced pairing + modern data); this table only serves as regression validation for the Srivastava extraction. Both conventions give p < 0.01, and the directional conclusion (systematic depolarization of cancer) is robust.

The pair-by-pair quality assessment is in the Level-1/2/3 provenance column and Notes 1–12 of Data 1.

---

## References (SI Data 2)

1. Srivastava P, Kane A, Harrison C, Levin M. A Meta-Analysis of Bioelectric Data in Cancer, Embryogenesis, and Regeneration. *Bioelectricity* 2021;3(1):42–67. PMID: 34476377
2. Binggeli R, Weinstein RC. Deficits in elevating membrane potential of rat fibrosarcoma cells after cell contact. *Cancer Res* 1985;45:235–241. PMID: 3965134
3. Marino AA, Iliev IG, Schwalke MA, et al. Association between cell membrane potential and breast cancer. *Tumour Biol* 1994;15:82–89. PMID: 8184256
4. Killion JJ. Electrical properties of normal and transformed mammalian cells. *Biophys J* 1984;45(3):523–528. PMID: 6324916

> **Note**: The complete sources and PMIDs of the remaining references cited in this table from Srivastava 2021 Supplementary Table S1 (Lash 1955, Tokuoka 1957, Jamakosmanovic 1968, Sand 1986, Picker 1981, Chowdhury 1973, Binggeli & Cameron 1980, Balitsky 1964, etc.) are in the Data 1 reference list #1–38.
