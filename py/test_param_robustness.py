# -*- coding: utf-8 -*-
"""Parameter robustness scan (3-equation convention): sweep 12 parameters to see the robustness of the main conclusions.

Main conclusions: P_max* (saddle-node), h_crit (reverse saddle-node = meridian open), low-health quiescence (h=0), bistability (h=1), high-health quiescence.
Expected (continuing from the I2 seven-equation version): low-health quiescence is structurally robust; high-health quiescence is fragile to "source too strong / clearance too weak";
    d_cns must = 0; γ has a narrow robust interval (≈0.35–0.4).
"""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import numpy as np

BASE = dict(u_star=0.4, k0=0.01, alpha_a=0.15, P_max=30.0, Z0=1.0, Zinf=0.05,
            phi=0.5, eps_u=5.0, h0=1.0, d0=1.0, gamma=0.5, d_cns=0.0,
            p0=2.0, e0=2.0, S0=0.1)

u_star = k0 = alpha_a = P_max = Z0 = Zinf = phi = eps_u = h0 = d0 = gamma = d_cns = p0 = e0 = S0 = None

def set_params(**kw):
    globals().update(kw)

set_params(**BASE)

def s_p(h): return h / (h + h0)
def Z(h): return (1.0 + (Zinf / Z0) * h) / (1.0 + h)
def k_tot(h): return k0 * phi * (1.0 + eps_u * h / (h + 1.0)) + k0 * (1.0 - phi) / Z(h)
def a_eff(h): return alpha_a * Z(h) / k_tot(h)
def P_eff(h, Pm=P_max): return Pm * s_p(h)

def d_of_S(S): return d_cns + (1.0 - d_cns) * gamma * S / (S + S0)
def p_of_d(d, h, Pm=P_max): return P_eff(h, Pm) * d / (d + d0)
def D_of_p(p, h):
    a = a_eff(h); f = p / (p + p0); e = max(0.0, a * p - u_star); return f * e / (e + e0)

def G(S, h, Pm=P_max): return D_of_p(p_of_d(d_of_S(S), h, Pm), h)

def _H(S, h, Pm): return G(S, h, Pm) - S

def fixed_points(h, Pm, n=3001):
    Sg = np.linspace(0.0, 1.0, n); H = np.array([_H(s, h, Pm) for s in Sg])
    roots = []
    if abs(H[0]) < 1e-12: roots.append(0.0)
    for i in range(n - 1):
        if H[i] * H[i + 1] < 0.0:
            lo, hi = Sg[i], Sg[i + 1]; Hlo = H[i]
            for _ in range(60):
                mid = 0.5 * (lo + hi); Hm = _H(mid, h, Pm)
                if Hlo * Hm <= 0.0: hi = mid
                else: lo = mid; Hlo = Hm
            roots.append(0.5 * (lo + hi))
    return sorted(set(round(r, 9) for r in roots))

def is_stable(S, h, Pm, eps=1e-5):
    if S < 1e-9: return True
    g = (G(S + eps, h, Pm) - G(S - eps, h, Pm)) / (2.0 * eps)
    return abs(g) < 1.0

def stable_nonzero(h, Pm):
    return [r for r in fixed_points(h, Pm) if r > 1e-9 and is_stable(r, h, Pm)]

def Pmax_star():
    for Pm in np.arange(0.0, 40.0, 0.5):
        if stable_nonzero(1.0, Pm): return Pm
    return None

def h_crit():
    for h in np.arange(60.0, 0.0, -0.5):
        if stable_nonzero(h, P_max): return h
    return None

def low_silent(): return len(stable_nonzero(0.0, P_max)) == 0
def bistable(): return len(stable_nonzero(1.0, P_max)) > 0

sweeps = {
    "d_cns": [0.0, 0.02, 0.05, 0.1],
    "gamma": [0.2, 0.3, 0.35, 0.4, 0.5, 0.6],
    "S0":    [0.02, 0.05, 0.1, 0.2, 0.5],
    "d0":    [0.3, 0.5, 1.0, 2.0, 4.0],
    "u_star":[0.2, 0.3, 0.4, 0.5, 0.6],
    "p0":    [0.5, 1.0, 2.0, 4.0, 8.0],
    "e0":    [0.5, 1.0, 2.0, 4.0, 8.0],
    "h0":    [0.2, 0.5, 1.0, 2.0, 5.0],
    "eps_u": [2.0, 3.0, 5.0, 8.0, 10.0],
    "phi":   [0.2, 0.35, 0.5, 0.65, 0.8],
    "Zinf":  [0.02, 0.03, 0.05, 0.08, 0.12],
    "P_max": [10, 20, 30, 40, 60],
}

def fmt(x):
    if isinstance(x, float):
        s = ("%.2f" % x).rstrip("0").rstrip(".")
        return s if s else "0"
    return str(x)

print("param     | value | P_max* | h_crit | low-h quiesc. | bistable | high-h quiesc.")
print("-" * 76)
for name, vals in sweeps.items():
    for v in vals:
        set_params(**BASE)
        set_params(**{name: v})
        pm = Pmax_star(); hc = h_crit()
        ls = low_silent(); bs = bistable()
        hs = hc is not None and hc < 60.0
        print("%-9s | %5s | %6s | %6s | %s | %s | %s" % (
            name, fmt(v),
            ("%.1f" % pm) if pm else "none",
            ("%.1f" % hc) if hc else "none",
            "✓" if ls else "✗", "✓" if bs else "✗", "✓" if hs else "✗"))
