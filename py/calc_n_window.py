# -*- coding: utf-8 -*-
"""Final criterion (h0=10, P_max=16, α_a=0.28125, γ=0.4): recheck the PSC window for n=0.5/1/2/3 (vectorized fast version).
has_psc ⟺ max_S (G(S)-S) > 0 (S=0 is always a fixed point, G(0)=0, G(1)<1, so a positive crossing implies a stable nonzero fixed point)."""
import sys
try: sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception: pass
import numpy as np

u_star=0.4; k0=0.01; alpha_a=0.28125; P_max=16.0; Z0=1.0; Zinf=0.05
phi=0.5; eps_u=5.0; h0=10.0; d0=1.0; d_cns=0.0; e0=2.0; S0=0.1; S_seed=0.8

def s_p(h): return h/(h+h0)
def Z(h): return (1.0+(Zinf/Z0)*h)/(1.0+h)
def k_tot(h): return k0*phi*(1.0+eps_u*h/(h+1.0))+k0*(1.0-phi)/Z(h)
def a_eff(h): return alpha_a*Z(h)/k_tot(h)
def P_eff(h): return P_max*s_p(h)

Sg = np.linspace(0.0, 1.0, 3001)
def H_of_S(h, gamma, n):
    d = gamma*Sg/(Sg+S0)                 # d_cns=0
    p = P_eff(h)*d/(d+d0)
    u = a_eff(h)*p
    x = np.maximum(0.0, u-u_star)
    G = x**n/(x**n+e0**n)
    return G - Sg

def has_psc(h, gamma, n):
    return H_of_S(h, gamma, n).max() > 1e-9

def G_at_seed(h, gamma, n):
    d = gamma*S_seed/(S_seed+S0)
    p = P_eff(h)*d/(d+d0)
    u = a_eff(h)*p
    x = max(0.0, u-u_star)
    return x**n/(x**n+e0**n)

def has_sense(h, gamma, n):
    return G_at_seed(h, gamma, n) > 1e-9

HMAX = 1e5
def bisect_flip(pred, lo, hi, want_high):
    """pred(lo)=want_high, pred(hi)=not want_high (or reversed) → return the flip point."""
    plo = pred(lo)
    for _ in range(90):
        mid = 0.5*(lo+hi)
        if pred(mid) == plo: lo = mid
        else: hi = mid
    return 0.5*(lo+hi)

def find_hcrit(gamma, n):
    if has_psc(HMAX, gamma, n):
        return None  # never quiescent
    hs = np.logspace(np.log10(HMAX), -3, 1200)
    for h in hs:
        if has_psc(h, gamma, n):
            return bisect_flip(lambda hh: has_psc(hh, gamma, n), h, HMAX, True)
    return None

def find_hlow(gamma, n):
    hs = np.logspace(-3, np.log10(HMAX), 1200)
    for h in hs:
        if has_psc(h, gamma, n):
            return bisect_flip(lambda hh: not has_psc(hh, gamma, n), 1e-9, h, False)
    return None

def find_hsense(gamma, n):
    hs = np.logspace(-3, np.log10(HMAX), 1200)
    for h in hs:
        if has_sense(h, gamma, n):
            return bisect_flip(lambda hh: not has_sense(hh, gamma, n), 1e-9, h, False)
    return None

Lam_inf = a_eff(1e6)*P_eff(1e6)/u_star
print("Λ(∞)=%.4f (Zinf=0.05 positive lower bound)"%Lam_inf)
print("="*80)
for gamma in [0.4, 0.8]:
    print("\nγ=%.1f:"%gamma)
    print("   %-6s %-12s %-12s %-14s %-12s"%("n","h_sense","h_low","h_crit","window width"))
    for n in [0.5,1.0,2.0,3.0]:
        hs=find_hsense(gamma,n); hl=find_hlow(gamma,n); hc=find_hcrit(gamma,n)
        w=(hc-hl) if (hl is not None and hc is not None) else None
        print("   %-6s %-12s %-12s %-14s %-12s"%(
            "%.1f"%n,
            "%.3f"%hs if hs is not None else "none",
            "%.3f"%hl if hl is not None else "none",
            "%.1f"%hc if hc is not None else "never quiescent",
            "%.1f"%w if w is not None else "—"))
