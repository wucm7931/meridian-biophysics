"""stat_rmp_cancer.py — recompute the RMP sign test of 04 main text §2.4.2 and Appendix 2, sub-appendices 2.3/2.4.

Data sources (v12/data/):
  - paired_rmp_data.csv              33 normal-cancer pairs (+1 HaCaT reference row, quality=Level-3 excluded)
  - srivastava_2021_extracted.csv    Srivastava 2021 extracted pairs (paired_id non-empty, GA2 inflammation control excluded)

Classification rule (row-by-row consistent with SI Table S1 "type" column):
  Let Δ = |RMP|_N − |RMP|_C (mV).
    depolarization (cancer R < normal R): Δ > 0        — no magnitude threshold on the depolarization side
    hyperpolarization (cancer R > normal R): Δ < −5 mV    — |Δ| > 5 mV and cancer more polarized
    tie: −5 ≤ Δ ≤ 0 mV                   — cancer slightly polarized but within noise (astrocyte I −1, LNCaP −3.5 mV)
  Sign test = two-sided exact binomial test (scipy.stats.binomtest, p0=0.5), tie pairs excluded.

Five p values in the main text (this script's output, iron rule ① digit-for-digit check):
  (1) main statistic       30/31 → 3.0e-8
  (2) Level-1 subset       13/13 → 2.4e-4
  (3) conservative merge   19/20 → 4.0e-5
  (4) Srivastava criterion A 26/28 → 3.0e-6
  (5) Srivastava criterion B 13/15 → 7.4e-3
"""

import os
from scipy.stats import binomtest

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, '..', 'data'))
PAIRED = os.path.join(DATA, 'paired_rmp_data.csv')
SRIV = os.path.join(DATA, 'srivastava_2021_extracted.csv')

# Conservative merge: pseudo-replicate groups sharing the same normal reference (cell_type exact match)
MERGE_GROUPS = {'胶质细胞', '肝细胞', '乳腺上皮', '黑色素细胞'}

# Valid paired IDs in the Srivastava extraction table (listed pair-by-pair in SI Table S2; GA2 is a gastritis inflammation control, excluded from the sign test)
VALID_PIDS = {'CX1', 'GA1', 'GA2', 'TH1', 'BR1', 'BR2', 'BR3', 'GL_all',
              'SK1', 'TH2', 'LI1', 'FI1', 'FI2', 'TH3', 'MU1', '3T1'}


def classify(delta_mv):
    """Δ = |RMP|_N − |RMP|_C (mV)."""
    if delta_mv > 0:
        return '去极化'
    if delta_mv < -5.0:
        return '超极化'
    return '持平'


def tally(pairs):
    """pairs: iterable of (normal_mV, cancer_mV) → (depolarization, hyperpolarization, tie)."""
    d = h = t = 0
    for n, c in pairs:
        k = classify(n - c)
        if k == '去极化':
            d += 1
        elif k == '超极化':
            h += 1
        else:
            t += 1
    return d, h, t


def sign_test(d, h, label):
    n = d + h
    p = binomtest(d, n, 0.5).pvalue
    print(f'{label}: depolarization {d}/{n} = {d/n:.1%}, hyperpolarization {h}/{n}, '
          f'non-tie n={n}, two-sided exact binomial p = {p:.2e}')
    return p


def load_paired():
    rows = []
    with open(PAIRED, encoding='utf-8') as f:
        for line in f:
            s = line.rstrip('\n').rstrip('\r')
            if not s or s.startswith('cell_type'):
                continue
            p = s.split(',')
            rows.append({
                'cell_type': p[0].strip(),
                'normal': float(p[1]),
                'cancer': float(p[2]),
                'quality': p[5].strip(),
            })
    return rows


def load_srivastava():
    groups = {}
    with open(SRIV, encoding='utf-8') as f:
        for line in f:
            s = line.rstrip('\n').rstrip('\r')
            if not s or s.startswith('species') or s.startswith('#'):
                continue
            p = s.split(',')
            pid = p[8].strip() if len(p) > 8 else ''
            if pid not in VALID_PIDS or pid == 'GA2':  # keep only valid paired rows (exclude rows without a paired reference)
                continue
            rtype = p[3].strip()
            rmp = float(p[4])
            g = groups.setdefault(pid, {'normal': None, 'cancers': []})
            if rtype == 'Normal':
                g['normal'] = rmp
            elif rtype == 'Cancer':
                g['cancers'].append(rmp)
    return groups


def main():
    print('=' * 72)
    print('RMP sign-test recomputation')
    print('=' * 72)

    # ---- 33-pair main statistic + Level-1 + conservative merge ----
    rows = load_paired()
    valid = [r for r in rows if r['quality'] != 'Level-3']  # exclude the HaCaT reference row
    assert len(valid) == 33, f'expected 33 valid pairs, got {len(valid)}'

    d, h, t = tally((r['normal'], r['cancer']) for r in valid)
    print(f'[main] 33 pairs: depolarization {d}, hyperpolarization {h}, tie {t}')
    assert (d, h, t) == (30, 1, 2), (d, h, t)
    sign_test(d, h, '(1) main statistic')

    l1 = [r for r in valid if r['quality'] == 'Level-1']
    d1, h1, t1 = tally((r['normal'], r['cancer']) for r in l1)
    print(f'[Level-1] {len(l1)} pairs: depolarization {d1}, hyperpolarization {h1}, tie {t1}')
    assert (len(l1), d1, h1, t1) == (14, 13, 0, 1), (len(l1), d1, h1, t1)
    sign_test(d1, h1, '(2) Level-1 subset')

    # Conservative merge: average groups sharing the same normal reference (glioma 5 grades, liver cancer 4 lines, breast 5 rows, melanoma 2 rows)
    merged = []
    grouped = {}
    for r in valid:
        if r['cell_type'] in MERGE_GROUPS:
            grouped.setdefault(r['cell_type'], []).append(r)
        else:
            merged.append(r)
    for ct, rs in grouped.items():
        merged.append({'normal': sum(x['normal'] for x in rs) / len(rs),
                       'cancer': sum(x['cancer'] for x in rs) / len(rs)})
    dm, hm, tm = tally((r['normal'], r['cancer']) for r in merged)
    print(f'[conservative merge] {len(merged)} pairs: depolarization {dm}, hyperpolarization {hm}, tie {tm}')
    assert (len(merged), dm, hm, tm) == (21, 19, 1, 1), (len(merged), dm, hm, tm)
    sign_test(dm, hm, '(3) conservative merge')

    # ---- Srivastava criterion A / B ----
    groups = load_srivastava()
    pairs_A = [(abs(g['normal']), abs(c)) for g in groups.values() for c in g['cancers']]
    dA, hA, tA = tally(pairs_A)
    print(f'[Srivastava criterion A row-wise] {len(pairs_A)} pairs: depolarization {dA}, hyperpolarization {hA}, tie {tA}')
    assert (len(pairs_A), dA, hA, tA) == (29, 26, 2, 1), (len(pairs_A), dA, hA, tA)
    sign_test(dA, hA, '(4) Srivastava criterion A')

    pairs_B = [(abs(g['normal']), abs(sum(g['cancers']) / len(g['cancers'])))
               for g in groups.values()]
    dB, hB, tB = tally(pairs_B)
    print(f'[Srivastava criterion B merged] {len(pairs_B)} pairs: depolarization {dB}, hyperpolarization {hB}, tie {tB}')
    assert (len(pairs_B), dB, hB, tB) == (15, 13, 2, 0), (len(pairs_B), dB, hB, tB)
    sign_test(dB, hB, '(5) Srivastava criterion B')

    # ---- Sensitivity analysis (SI Table S2: exclude the 2 contested hyperpolarized pairs TH1, 3T1) ----
    sens_A = [(n, c) for (n, c) in pairs_A if classify(n - c) != '超极化']
    dSA, hSA, tSA = tally(sens_A)
    assert (dSA, hSA, tSA) == (26, 0, 1), (dSA, hSA, tSA)
    sign_test(dSA, hSA, '(sensitivity) criterion A excluding hyperpolarized pairs')

    sens_B = [(n, c) for (n, c) in pairs_B if classify(n - c) != '超极化']
    dSB, hSB, tSB = tally(sens_B)
    assert (dSB, hSB, tSB) == (13, 0, 0), (dSB, hSB, tSB)
    sign_test(dSB, hSB, '(sensitivity) criterion B excluding hyperpolarized pairs')

    print('=' * 72)
    print('All counts match the main text (assertions passed).')


if __name__ == '__main__':
    main()
