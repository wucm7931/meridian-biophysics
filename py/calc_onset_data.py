# -*- coding: utf-8 -*-
"""Systematic comparison of high/low PSC onset (h_low): sweep n × γ × h0 (+ Zinf fixed = 0.05).
Outputs per combination: no-sensation→sensation onset h_sense, PSC onset h_low, PSC end h_crit, window width.
Purpose: find which knob decides "high/low PSC onset"."""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import numpy as np

u_star=0.4; k0=0.01; alpha_a=0.28125; P_max=16.0; Z0=1.0; Zinf=0.05  # P_max=16 aligns with point ⑤ P_eff=16; α_a=0.15*30/16, α_a·P_max=4.5 unchanged→boundary unchanged
phi=0.5; eps_u=5.0; d0=1.0; d_cns=0.0; e0=2.0; S0=0.1; S_seed=0.8

def make(h0):
    def s_p(h): return h/(h+h0)
    def Z(h): return (1.0+(Zinf/Z0)*h)/(1.0+h)
    def k_tot(h): return k0*phi*(1.0+eps_u*h/(h+1.0))+k0*(1.0-phi)/Z(h)
    def a_eff(h): return alpha_a*Z(h)/k_tot(h)
    def P_eff(h): return P_max*s_p(h)
    def d_of_S(S,gamma): return d_cns+(1.0-d_cns)*gamma*S/(S+S0)
    def p_of_d(d,h): return P_eff(h)*d/(d+d0)
    def D_of_u(u,n):
        x=max(0.0,u-u_star); return x**n/(x**n+e0**n)
    def G(S,h,gamma,n): return D_of_u(a_eff(h)*p_of_d(d_of_S(S,gamma),h), n)
    def _H(S,h,gamma,n): return G(S,h,gamma,n)-S
    return locals()

def analyze(h0, gamma, n, hmax=400.0):
    E=make(h0)
    # log + linear mixed grid, fine sweep of low health
    hs=sorted(set(
        [0.0,0.001,0.002,0.003,0.005,0.008,0.01,0.015,0.02,0.03,0.05,0.08,0.1,0.15,0.2,0.3,0.5,0.8,
         1.0,2.0,5.0,10.0,20.0,50.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0]
    ))
    hs=[h for h in hs if h<=hmax]
    def fps(h):
        Sg=np.linspace(0,1,3001); H=np.array([E['_H'](s,h,gamma,n) for s in Sg]); roots=[]
        for i in range(len(Sg)-1):
            if H[i]*H[i+1]<0:
                lo,hi=Sg[i],Sg[i+1]; Hlo=H[i]
                for _ in range(50):
                    mid=0.5*(lo+hi); Hm=E['_H'](mid,h,gamma,n)
                    if Hlo*Hm<=0: hi=mid
                    else: lo=mid; Hlo=Hm
                roots.append(0.5*(lo+hi))
        return sorted(set(round(r,8) for r in roots))
    def is_stable(S,h,eps=1e-5):
        if S<1e-9: return True
        g=(E['G'](S+eps,h,gamma,n)-E['G'](S-eps,h,gamma,n))/(2*eps)
        return abs(g)<1.0
    def has_psc(h): return any(is_stable(r,h) for r in fps(h) if r>1e-9)
    def has_sense(h): return E['G'](S_seed,h,gamma,n)>1e-9
    # sequential sweep, record state flip points
    h_sense=None; h_low=None; h_crit=None
    prev_sense=None; prev_psc=None
    for h in hs:
        s=has_sense(h); p=has_psc(h)
        if prev_sense is None:
            prev_sense=s; prev_psc=p
            continue
        if not prev_sense and s and h_sense is None: h_sense=h
        if not prev_psc and p and h_low is None: h_low=h
        if prev_psc and not p and h_crit is None: h_crit=h
        prev_sense=s; prev_psc=p
    return h_sense, h_low, h_crit

print("Zinf=0.05 fixed. Each cell = 'PSC onset h_low / end h_crit / window width', also marks sensation onset h_sense.")
print("h0 = source-health half-saturation (↑→P_eff rises slowly); γ = retreat depth; n = Hill exponent.")
print("="*100)

for n in [0.5,1.0,2.0,3.0]:
    print("\n■ n=%.1f"%n)
    print("  %-20s %-12s %-12s %-12s %-12s"%("h0, γ","h_sense","h_low(onset)","h_crit(end)","window width"))
    for h0 in [1.0,10.0,30.0]:
        for g in [0.3,0.5,0.8]:
            hs,hl,hc=analyze(h0,g,n)
            w = (hc-hl) if (hl is not None and hc is not None) else None
            print("  h0=%-5.1f γ=%.1f   %-12s %-12s %-12s %-12s"%(
                h0,g,
                "%.3f"%hs if hs else "none",
                "%.3f"%hl if hl else "none (no PSC at all)",
                "%.1f"%hc if hc else "none",
                "%.1f"%w if w else "—"))
