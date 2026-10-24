"""Fig 2.4-1: configurational entropy & |E_K| vs K+ intracellular fraction f_i,
comparing extracellular volume fractions alpha = 0.5 (equal volume, shown first)
and alpha = 0.2 (physiological).

Uses the paper's unequal-volume two-compartment forms:
    S_config/k_B = -[ f_i ln(f_i/(1-a)) + f_o ln(f_o/a) ]          (max 0 at f_i = 1-a)
    E_K = (RT/F) ln( (f_o/f_i) * ((1-a)/a) )                        (zero at f_i = 1-a)

Original D:\\9-median_v2\\psc_model\\fig_entropy_rmp.py(α-corrected version), 2026-10-20 copied into
v12/py/, only changed the output directory to v12/images/ (not a placeholder figure; Fig 2.4-1 is a script-computed figure, do not hand-draw).
"""
import os
import matplotlib.pyplot as plt
import numpy as np

RT_F = 26.7  # mV at 37 C

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7.5), sharex=True)

f_i = np.linspace(0.01, 0.99, 600)
f_o = 1.0 - f_i

# (alpha, color) — alpha = 0.5 drawn and labelled first
alphas = [(0.5, '#e67e22'), (0.2, '#2a6496')]

for a, color in alphas:
    lbl = r'$\alpha$ = %.1f' % a
    # configurational entropy per ion (units of k_B), referenced to equal concentration
    S = -(f_i * np.log(f_i / (1.0 - a)) + f_o * np.log(f_o / a))
    ax1.plot(f_i, S, color=color, lw=2.5, label=lbl)

    # concentration ratio [K+]o/[K+]i  -> Nernst |E_K|
    r = (f_o / f_i) * ((1.0 - a) / a)
    ax2.plot(f_i, np.abs(RT_F * np.log(r)), color=color, lw=2.5, label=lbl)

# ---- top panel: equal-concentration (max entropy) reference points ----
ax1.axhline(0.0, color='black', lw=0.8, alpha=0.35)
for a, color in alphas:
    feq = 1.0 - a
    ax1.plot(feq, 0.0, 'o', color=color, ms=8, zorder=6)
    ax1.axvline(feq, color=color, ls=':', lw=1, alpha=0.55)
ax1.set_ylabel(r'$S_{config}\,/\,k_B$', fontsize=12)
ax1.set_ylim(-1.8, 0.25)
ax1.legend(fontsize=10, loc='lower right', frameon=False)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# ---- bottom panel ----
def rmp_at(fi, a):
    r = ((1.0 - fi) / fi) * ((1.0 - a) / a)
    return abs(RT_F * np.log(r))

for a, color in alphas:
    feq = 1.0 - a
    ax2.plot(feq, 0.0, 'o', color=color, ms=8, zorder=6)
    ax2.axvline(feq, color=color, ls=':', lw=1, alpha=0.55)

# physiological states on the alpha = 0.2 curve
states = [(0.8, 'Equal\nconcentration'),
          (0.947, 'Depolarized\ncells\n(~40 mV)'),
          (0.981, 'Mature\nneuron\n(~69 mV)')]
for fi, label in states:
    v = rmp_at(fi, 0.2)
    ax2.plot(fi, v, 'o', color='#2a6496', ms=8, zorder=6)
    ax2.annotate(label, xy=(fi, v), fontsize=9, ha='center', va='bottom',
                 xytext=(0, 10), textcoords='offset points',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.85))

ax2.set_xlabel(r'$f_i$ (K$^+$ intracellular fraction)', fontsize=12)
ax2.set_ylabel(r'$|E_K|$ = $|$RMP$|$ (mV)', fontsize=12)
ax2.set_ylim(-5, 110)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

plt.tight_layout(pad=1.5)

here = os.path.dirname(os.path.abspath(__file__))
out_dir = os.path.normpath(os.path.join(here, '..', 'images'))
os.makedirs(out_dir, exist_ok=True)
out = os.path.join(out_dir, 'fig03_entropy.png')
plt.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
print('Saved:', out)
