"""Figure 2.5-1: Three-zone Hormesis dose-response curve (English version).

Three zones (reconstructed concept):
  Zone F (under-stimulation): below MES -> alpha2 inhibition -> R passively declines -> A-
  Zone E (moderate): MES to damage threshold -> alpha2 upregulation -> E-ratchet up -> A+
  Excess E (over-stress): above damage threshold -> adaptation reversal -> R falls back A- (abnormal, reversible); more extreme -> C/D

Inverted-U curve (rise -> plateau -> decline); both ends A-, middle A+.

Originally D:/9-median_v2/psc_model/fig_three_zone_hormesis.py (Chinese), copied to v12/py/ 2026-10-20;
2026-10-21 translated to English + recolored (blue/green/red). Qualitative schematic; MES/DAMAGE are hand-set illustrative values.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def setup_style():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'axes.unicode_minus': False,
        'font.size': 11,
    })


def main():
    setup_style()
    fig, ax = plt.subplots(figsize=(10.5, 5.8))

    MES = 1.0       # minimal effective stimulus (Zone F / Zone E boundary)
    DAMAGE = 3.5    # damage threshold (Zone E / Excess E boundary)
    XMAX = 5.5

    x = np.linspace(0, XMAX, 600)
    y = np.empty_like(x)
    m_f = x <= MES
    m_e = (x > MES) & (x <= DAMAGE)
    m_h = x > DAMAGE
    # F: low stimulus, slow rise
    y[m_f] = 0.12 + 0.28 * (x[m_f] / MES)
    # E: adaptation rising to plateau
    y[m_e] = 0.40 + 0.55 * (x[m_e] - MES) / (DAMAGE - MES)
    # Excess E: over-threshold decline (not to 0: abnormal, reversible; only more extreme -> C/D)
    y[m_h] = 0.95 - 0.70 * (x[m_h] - DAMAGE) / (XMAX - DAMAGE)
    y = np.clip(y, 0, 1)

    C_F, C_E, C_H = '#4C72B0', '#55A868', '#FF8000'

    # three colored bands
    ax.fill_between(x[m_f], 0, y[m_f], alpha=0.20, color=C_F)
    ax.fill_between(x[m_e], 0, y[m_e], alpha=0.20, color=C_E)
    ax.fill_between(x[m_h], 0, y[m_h], alpha=0.20, color=C_H)

    # curve
    ax.plot(x, y, 'k-', lw=2.6, zorder=5)

    # threshold dashed lines
    ax.axvline(MES, color='#888888', ls=':', lw=1.4, zorder=4)
    ax.axvline(DAMAGE, color='#888888', ls=':', lw=1.4, zorder=4)

    # threshold labels (top)
    ax.text(MES, 1.09, 'MES (minimal effective stimulus)', ha='center', va='bottom',
            fontsize=10.5, color='#333333')
    ax.text(DAMAGE, 1.09, 'Damage threshold', ha='center', va='bottom',
            fontsize=10.5, color='#333333')

    # Zone F label
    ax.text(MES / 2, 0.08, 'Zone F', ha='center', va='bottom',
            fontsize=14, fontweight='bold', color=C_F)
    ax.text(MES / 2, 0.04, 'Under-stimulation \u2192 A\u2212', ha='center', va='bottom',
            fontsize=9.5, color=C_F)

    # Zone E label
    ax.text((MES + DAMAGE) / 2, 0.32, 'Zone E', ha='center', va='bottom',
            fontsize=14, fontweight='bold', color=C_E)
    ax.text((MES + DAMAGE) / 2, 0.24, 'Moderate stress \u2192 E-ratchet \u2192 A+',
            ha='center', va='bottom', fontsize=9.5, color=C_E)

    # Excess E label
    cx = (DAMAGE + XMAX) / 2
    ax.text(cx, 0.32, 'Excess E', ha='center', va='bottom',
            fontsize=14, fontweight='bold', color=C_H)
    ax.text(cx, 0.24, 'Excessive stress \u2192 A\u2212', ha='center', va='bottom',
            fontsize=9.5, color=C_H)
    ax.text(cx, 0.17, '(abnormal, reversible;', ha='center', va='bottom',
            fontsize=9.5, color=C_H)
    ax.text(cx, 0.10, 'more extreme \u2192 C/D)', ha='center', va='bottom',
            fontsize=9.5, color=C_H)

    # axes
    ax.set_xlim(0, XMAX)
    ax.set_ylim(0, 1.22)
    ax.set_xlabel('Stimulus intensity (dose)', fontsize=12.5, labelpad=8)
    ax.set_ylabel('Tissue function level', fontsize=12.5, labelpad=8)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()
    out = Path(__file__).resolve().parent.parent / 'images' / 'fig05_hormesis.png'
    fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    print(f'Saved: {out}')


if __name__ == '__main__':
    main()
