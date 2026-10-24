# -*- coding: utf-8 -*-
"""Decay-segment count vs Λ: confirm the relay is determined solely by the single parameter Λ = a·P_eff/u*, with furthest distance = f(Λ*-Λ).

⚠️ Current values come from Part 1 and Part 2 of this script — Λ*(n=2) and the "gap→segments" table (criterion-independent);
   Part 3 (health-layer comparison) uses the P_max=30 old criterion and has been superseded by `calc_n_window.py` (P_max=16); do not take values from it."""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

u_star=0.4; gamma=0.4; S0=0.1; d0=1.0; e0=2.0; n=2.0; S_seed=0.8

def d_of_S(S): return gamma*S/(S+S0)
def u_of_S(S, Lam): return Lam*u_star * d_of_S(S)/(d_of_S(S)+d0)
def D_of_u(u):
    x=max(0.0,u-u_star); return x**n/(x**n+e0**n)
def G(S, Lam): return D_of_u(u_of_S(S, Lam))

# ---- 1) find the saddle-node Λ*: where a stable nonzero fixed point first appears ----
def has_stable_nonzero(Lam):
    # scan the S grid for stable zeros of G(S)-S
    Sg=[i/2000.0 for i in range(2001)]
    H=[G(s,Lam)-s for s in Sg]
    roots=[]
    for i in range(len(Sg)-1):
        if H[i]*H[i+1]<0:
            lo,hi=Sg[i],Sg[i+1]; Hlo=H[i]
            for _ in range(60):
                mid=0.5*(lo+hi); Hm=G(mid,Lam)-mid
                if Hlo*Hm<=0: hi=mid
                else: lo=mid; Hlo=Hm
            roots.append(0.5*(lo+hi))
    for r in roots:
        if r>1e-9:
            eps=1e-5; g=(G(r+eps,Lam)-G(r-eps,Lam))/(2*eps)
            if abs(g)<1.0: return True
    return False

lo,hi=1.0,40.0
for _ in range(80):
    mid=0.5*(lo+hi)
    if has_stable_nonzero(mid): hi=mid
    else: lo=mid
Lam_star=hi
print("saddle-node Λ*(n=2,γ=0.4,S0=0.1,d0=1.0,u*=0.4,e0=2.0) = %.4f"%Lam_star)
print()

# ---- 2) decay-segment count vs Λ ----
def decay_segments(Lam):
    S=S_seed
    for k in range(10000):
        u=u_of_S(S,Lam)
        if u < u_star:  # this domain has fallen below threshold -> next domain S=0 quenches
            return k   # k = number of domains sensed (including the seed domain)
        S=G(S,Lam)
    return 10000  # self-sustaining

print("seed S=0.8. decay-segment count = f(Λ), Λ swept from far below Λ* up to Λ*:")
print("%-8s %-8s %-10s"%("Λ","Λ*-Λ","segments(≈cm)"))
print("-"*30)
for Lam in [Lam_star-8, Lam_star-6, Lam_star-4, Lam_star-2, Lam_star-1,
            Lam_star-0.5, Lam_star-0.2, Lam_star-0.1, Lam_star-0.05, Lam_star-0.01]:
    k=decay_segments(Lam)
    print("%-8.2f %-8.2f %-10s"%(Lam, Lam_star-Lam, k if k<10000 else "self-sustaining(∞)"))

# ---- 3) comparison: Λ and segment count for the real h, checking whether the low/high health branches lie on the same curve ----
print()
print("comparison (health functions locked): low-health branch (Λ rising) and high-health branch (Λ falling) with equal segment count = equal Λ.")
u_star2=0.4; k0=0.01; alpha_a=0.15; P_max=30.0; Z0=1.0; Zinf=0.05
phi=0.5; eps_u=5.0; h0=10.0
def s_p(h): return h/(h+h0)
def Z(h): return (1.0+(Zinf/Z0)*h)/(1.0+h)
def k_tot(h): return k0*phi*(1.0+eps_u*h/(h+1.0))+k0*(1.0-phi)/Z(h)
def a_eff(h): return alpha_a*Z(h)/k_tot(h)
def P_eff(h): return P_max*s_p(h)
def Lam_h(h): return a_eff(h)*P_eff(h)/u_star2

print("%-8s %-8s %-10s"%("h","Λ(h)","segments"))
for h in [0.10,0.20,0.30,0.35,0.40, 5.0,8.0,20.0,50.0]:
    L=Lam_h(h)
    print("%-8.2f %-8.2f %-10s"%(h, L, decay_segments(L) if decay_segments(L)<10000 else "self-sustaining"))
