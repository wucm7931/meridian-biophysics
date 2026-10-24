# meridian-biophysics

Companion code and data for:

> **From propagated sensation along the meridians to a quantifiable traditional Chinese medicine: a unified biophysical model**

Preprint: _bioRxiv — link to be added_

## Overview

Propagated sensation along the meridians (PSC) is a reproducible physiological phenomenon, but the meridians have lacked a unified anatomical basis and a quantitative theory. The paper proposes that PSC arises from interstitial K⁺ accumulation driven by sympathetic withdrawal and feedback-based parallel charging along sympathetic nerves, arterioles, and sensory terminals. Centered on normalized membrane-potential dynamics, it builds a six-mode framework (resting, excitation, proliferation, apoptosis, and ascending/descending adaptation) that gives measurable electrophysiological definitions of qi, yin–yang, and supporting-the-upright/dispelling-pathogen, and that explains PSC velocity, lesion tropism, and state transitions in cancer and aging.

The manuscript text is hosted with the preprint/journal; this repository holds the reproducible code, data, and figures.

## Repository contents

- `py/` — model, figure, and statistical scripts (one script per computation). See `py/README.md` for the "derived number → script" map.
- `data/` — supplementary tables (`SI_Table_S1–S4.md`) and the raw sign-test datasets (`paired_rmp_data.csv`, `srivastava_2021_extracted.csv`) used by `py/stat_rmp_cancer.py`.
- `images/` — figures 1–6 (300-DPI PNG). Figures 3 and 5 are script-computed and regenerable via `py/fig_entropy_rmp.py` and `py/fig_three_zone_hormesis.py`.
- `refs/` — bibliography: human-readable `final_references.md` and machine-readable `final_references.tsv`.

## Requirements and usage

- Python 3.8+ with `numpy`, `scipy`, `matplotlib`.
- Scripts are standalone; run them from the repository root, e.g. `python py/model_source_a.py` or `python py/stat_rmp_cancer.py`.

## License

- Code (`py/`) — MIT License (see `LICENSE`).
- Data, figures, and text — [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
