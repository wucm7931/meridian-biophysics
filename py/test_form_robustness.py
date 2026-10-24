# -*- coding: utf-8 -*-
"""Functional-form robustness scan (3-equation convention): how swapping the forms of
d(S), p(d), D(p) affects the main conclusions.

d(S) ∈ {hill, linear, step}; p(d) ∈ {hill, linear}; D(p) ∈ {thresh_hill, thresh_exp, sigmoid}
18 combinations in total, printing the main conclusion for each combination:
    P_max* (saddle-node), h_crit (reverse saddle-node = meridian open), low-health quiescence (h=0), bistability (h=1), high-health quiescence.
Expected: d/p forms are robust as long as they are monotonic + saturating; the hard threshold of D(p) is load-bearing (sigmoid without hard threshold → loses both quiescence and bistability).
"""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import numpy as np

# parameters (same as fixed_point_3eq.py)
u_star = 0.4; k0 = 0.01; alpha_a = 0.15; P_max = 30.0; Z0 = 1.0; Zinf = 0.05
phi = 0.5; eps_u = 5.0; h0 = 1.0; d0 = 1.0; gamma = 0.5; d_cns = 0.0
p0 = 2.0; e0 = 2.0; S0 = 0.1

def s_p(h): return h / (h + h0)
def Z(h): return (1.0 + (Zinf / Z0) * h) / (1.0 + h)
def k_tot(h): return k0 * phi * (1.0 + eps_u * h / (h + 1.0)) + k0 * (1.0 - phi) / Z(h)
def a_eff(h): return alpha_a * Z(h) / k_tot(h)
def P_eff(h, Pm=P_max): return Pm * s_p(h)

def d_of_S(S, form):
    if form == "hill": return d_cns + (1.0 - d_cns) * gamma * S / (S + S0)
    if form == "linear": return min(1.0, d_cns + (1.0 - d_cns) * gamma * S)
    if form == "step": return d_cns if S < S0 else d_cns + (1.0 - d_cns) * gamma
    raise ValueError(form)

def p_of_d(d, h, Pm=P_max, form="hill"):
    Pe = P_eff(h, Pm)
    if form == "hill": return Pe * d / (d + d0)
    if form == "linear": return Pe * min(1.0, d)
    raise ValueError(form)

def D_of_p(p, h, form):
    a = a_eff(h)
    if form == "thresh_hill":
        f = p / (p + p0); e = max(0.0, a * p - u_star); return f * e / (e + e0)
    if form == "thresh_exp":
        e = max(0.0, a * p - u_star); return (1.0 - np.exp(-p / p0)) * e / (e + e0)
    if form == "sigmoid":
        return 1.0 / (1.0 + np.exp(-(a * p - u_star) / e0))
    raise ValueError(form)

def G(S, h, Pm, forms):
    return D_of_p(p_of_d(d_of_S(S, forms[0]), h, Pm, forms[1]), h, forms[2])

def _H(S, h, Pm, forms): return G(S, h, Pm, forms) - S

def fixed_points(h, Pm, forms, n=3001):
    Sg = np.linspace(0.0, 1.0, n); H = np.array([_H(s, h, Pm, forms) for s in Sg])
    roots = []
    if abs(H[0]) < 1e-12: roots.append(0.0)
    for i in range(n - 1):
        if H[i] * H[i + 1] < 0.0:
            lo, hi = Sg[i], Sg[i + 1]; Hlo = H[i]
            for _ in range(60):
                mid = 0.5 * (lo + hi); Hm = _H(mid, h, Pm, forms)
                if Hlo * Hm <= 0.0: hi = mid
                else: lo = mid; Hlo = Hm
            roots.append(0.5 * (lo + hi))
    return sorted(set(round(r, 9) for r in roots))

def is_stable(S, h, Pm, forms, eps=1e-5):
    if S < 1e-9: return True
    g = (G(S + eps, h, Pm, forms) - G(S - eps, h, Pm, forms)) / (2.0 * eps)
    return abs(g) < 1.0

def stable_nonzero(h, Pm, forms):
    return [r for r in fixed_points(h, Pm, forms) if r > 1e-9 and is_stable(r, h, Pm, forms)]

def Pmax_star(forms):
    for Pm in np.arange(0.0, 40.0, 0.2):
        if stable_nonzero(1.0, Pm, forms): return Pm
    return None

def h_crit(forms):
    for h in np.arange(60.0, 0.0, -0.2):
        if stable_nonzero(h, P_max, forms): return h
    return None

d_forms = ["hill", "linear", "step"]
p_forms = ["hill", "linear"]
D_forms = ["thresh_hill", "thresh_exp", "sigmoid"]

print("forms (d, p, D)          | P_max* | h_crit | low-h quiesc. | bistable | high-h quiesc.")
print("-" * 78)
for df in d_forms:
    for pf in p_forms:
        for Df in D_forms:
            forms = (df, pf, Df)
            pm = Pmax_star(forms); hc = h_crit(forms)
            low_silent = len(stable_nonzero(0.0, P_max, forms)) == 0
            bistable = len(stable_nonzero(1.0, P_max, forms)) > 0 and abs(G(0, 1.0, P_max, forms)) < 1e-9
            high_silent = hc is not None and hc < 60.0
            print("%-26s | %6s | %6s | %s | %s | %s" % (
                "(%s,%s,%s)" % (df, pf, Df),
                ("%.1f" % pm) if pm else "none",
                ("%.1f" % hc) if hc else "none",
                "✓" if low_silent else "✗",
                "✓" if bistable else "✗",
                "✓" if high_silent else "✗"))
