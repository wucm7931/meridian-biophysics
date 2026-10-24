# -*- coding: utf-8 -*-
"""Test the power forms of S(p)=D(p): whether the main conclusion is robust.

Forms (the p→S relation of D, with x=a·p−u*, y=a·p):
  pow1      : x_+/(x_+ + e0)                             (thresholded Hill n=1)
  pow2      : x_+^2/(x_+^2 + e0^2)
  pow3      : x_+^3/(x_+^3 + e0^3)
  pow_half  : x_+^0.5/(x_+^0.5 + e0^0.5)
  pure_pow1 : y/(y + e0)                                  (pure power n=1, no hard threshold)
  pure_pow2 : y^2/(y^2 + e0^2)                            (pure power n=2, no hard threshold)
"""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import numpy as np

u_star = 0.4; k0 = 0.01; alpha_a = 0.15; P_max = 30.0; Z0 = 1.0; Zinf = 0.05
phi = 0.5; eps_u = 5.0; h0 = 1.0; d0 = 1.0; gamma = 0.5; d_cns = 0.0
p0 = 2.0; e0 = 2.0; S0 = 0.1

def s_p(h): return h / (h + h0)
def Z(h): return (1.0 + (Zinf / Z0) * h) / (1.0 + h)
def k_tot(h): return k0 * phi * (1.0 + eps_u * h / (h + 1.0)) + k0 * (1.0 - phi) / Z(h)
def a_eff(h): return alpha_a * Z(h) / k_tot(h)
def P_eff(h, Pm=P_max): return Pm * s_p(h)

def d_of_S(S): return d_cns + (1.0 - d_cns) * gamma * S / (S + S0)
def p_of_d(d, h, Pm=P_max): return P_eff(h, Pm) * d / (d + d0)

def D_of_p(p, h, form):
    a = a_eff(h)
    x = a * p - u_star
    y = a * p
    if form == "pow1":
        x = max(0.0, x); return x / (x + e0)
    if form == "pow2":
        x = max(0.0, x); return x * x / (x * x + e0 * e0)
    if form == "pow3":
        x = max(0.0, x); return x ** 3 / (x ** 3 + e0 ** 3)
    if form == "pow_half":
        x = max(0.0, x); s = x ** 0.5; return s / (s + e0 ** 0.5)
    if form == "pure_pow1":
        y = max(0.0, y); return y / (y + e0)
    if form == "pure_pow2":
        y = max(0.0, y); return y * y / (y * y + e0 * e0)
    if form == "log_thresh":   # thresholded log (Weber–Fechner type): hard threshold + slow log increase, hard cutoff at x=e0
        x = max(0.0, (a * p - u_star) / e0)
        return min(1.0, np.log(1.0 + x) / np.log(2.0))
    if form == "log_pure":     # thresholdless log: no hard threshold
        x = max(0.0, a * p / e0)
        return min(1.0, np.log(1.0 + x) / np.log(2.0))
    raise ValueError(form)

def G(S, h, Pm, form): return D_of_p(p_of_d(d_of_S(S), h, Pm), h, form)
def _H(S, h, Pm, form): return G(S, h, Pm, form) - S

def fixed_points(h, Pm, form, n=3001):
    Sg = np.linspace(0.0, 1.0, n); H = np.array([_H(s, h, Pm, form) for s in Sg])
    roots = []
    if abs(H[0]) < 1e-12: roots.append(0.0)
    for i in range(n - 1):
        if H[i] * H[i + 1] < 0.0:
            lo, hi = Sg[i], Sg[i + 1]; Hlo = H[i]
            for _ in range(60):
                mid = 0.5 * (lo + hi); Hm = _H(mid, h, Pm, form)
                if Hlo * Hm <= 0.0: hi = mid
                else: lo = mid; Hlo = Hm
            roots.append(0.5 * (lo + hi))
    return sorted(set(round(r, 9) for r in roots))

def is_stable(S, h, Pm, form, eps=1e-5):
    if S < 1e-9: return True
    g = (G(S + eps, h, Pm, form) - G(S - eps, h, Pm, form)) / (2.0 * eps)
    return abs(g) < 1.0

def stable_nonzero(h, Pm, form):
    return [r for r in fixed_points(h, Pm, form) if r > 1e-9 and is_stable(r, h, Pm, form)]

def Pmax_star(form):
    for Pm in np.arange(0.0, 40.0, 0.2):
        if stable_nonzero(1.0, Pm, form): return Pm
    return None

def h_crit(form):
    for h in np.arange(60.0, 0.0, -0.5):
        if stable_nonzero(h, P_max, form): return h
    return None

forms = ["pow1", "pow2", "pow3", "pow_half", "pure_pow1", "pure_pow2", "log_thresh", "log_pure"]
print("D(p) form           | P_max* | h_crit | low quiesc. | bistable | high quiesc.")
print("-" * 66)
for f in forms:
    pm = Pmax_star(f); hc = h_crit(f)
    ls = len(stable_nonzero(0.0, P_max, f)) == 0
    bs = len(stable_nonzero(1.0, P_max, f)) > 0
    hs = hc is not None and hc < 60.0
    print("%-18s | %6s | %6s | %s | %s | %s" % (
        f, ("%.1f" % pm) if pm else "none", ("%.1f" % hc) if hc else "none",
        "✓" if ls else "✗", "✓" if bs else "✗", "✓" if hs else "✗"))
