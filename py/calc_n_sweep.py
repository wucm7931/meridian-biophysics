# -*- coding: utf-8 -*-
"""n comparison (γ=0.4 unified criterion): Λ*, S*(Λ=20), h_low/h_crit/window, short-range segment-count curves.

⚠️ Current values come from Part 1 (mapping layer) of this script — the Λ* and S*(Λ=20) columns (criterion-independent);
   Part 2 (health layer h_low/h_crit) uses the P_max=30 old criterion and has been superseded by `calc_n_window.py`
   (P_max=16 final criterion); do not take values from it."""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

u_star=0.4; gamma=0.4; S0=0.1; d0=1.0; e0=2.0; S_seed=0.8

def d_of_S(S): return gamma*S/(S+S0)
def u_of_S(S, Lam): return Lam*u_star * d_of_S(S)/(d_of_S(S)+d0)
def D_of_u(u,n):
    x=max(0.0,u-u_star); return x**n/(x**n+e0**n)
def G(S, Lam, n): return D_of_u(u_of_S(S,Lam), n)

def has_stable_nonzero(Lam, n):
    Sg=[i/2000.0 for i in range(2001)]
    H=[G(s,Lam,n)-s for s in Sg]; roots=[]
    for i in range(len(Sg)-1):
        if H[i]*H[i+1]<0:
            lo,hi=Sg[i],Sg[i+1]; Hlo=H[i]
            for _ in range(60):
                mid=0.5*(lo+hi); Hm=G(mid,Lam,n)-mid
                if Hlo*Hm<=0: hi=mid
                else: lo=mid; Hlo=Hm
            roots.append(0.5*(lo+hi))
    for r in roots:
        if r>1e-9:
            eps=1e-5; g=(G(r+eps,Lam,n)-G(r-eps,Lam,n))/(2*eps)
            if abs(g)<1.0: return True
    return False

def bisect_Lamstar(n):
    lo,hi=1.0,60.0
    for _ in range(80):
        mid=0.5*(lo+hi)
        if has_stable_nonzero(mid,n): hi=mid
        else: lo=mid
    return hi

def S_star_at(Lam,n):
    Sg=[i/4000.0 for i in range(4001)]
    H=[G(s,Lam,n)-s for s in Sg]; roots=[]
    for i in range(len(Sg)-1):
        if H[i]*H[i+1]<0:
            lo,hi=Sg[i],Sg[i+1]; Hlo=H[i]
            for _ in range(60):
                mid=0.5*(lo+hi); Hm=G(mid,Lam,n)-mid
                if Hlo*Hm<=0: hi=mid
                else: lo=mid; Hlo=Hm
            roots.append(0.5*(lo+hi))
    # take the largest stable nonzero root
    best=0.0
    for r in roots:
        if r>1e-9:
            eps=1e-5; g=(G(r+eps,Lam,n)-G(r-eps,Lam,n))/(2*eps)
            if abs(g)<1.0 and r>best: best=r
    return best

def decay_segments(Lam,n):
    S=S_seed
    for k in range(10000):
        if u_of_S(S,Lam) < u_star: return k
        S=G(S,Lam,n)
    return 10000

print("=== 1) mapping layer (γ=0.4,S0=0.1,d0=1.0,u*=0.4,e0=2.0) ===")
print("%-6s %-10s %-12s %-14s %-14s"%("n","Λ*","S*(Λ=20)","Λ*-1 gap segments","Λ*-0.05 gap segments"))
for n in [0.5,1.0,2.0,3.0]:
    Ls=bisect_Lamstar(n)
    Ss=S_star_at(20.0,n)
    k1=decay_segments(Ls-1.0,n)
    k2=decay_segments(Ls-0.05,n)
    print("%-6.1f %-10.3f %-12.3f %-14s %-14s"%(n,Ls,Ss,k1,k2))

print()
print("=== 2) health layer (locked h0=10,γ=0.4,Zinf=0.05; P_max/α_a placeholder) ===")
u_star2=0.4; k0=0.01; alpha_a=0.15; P_max=30.0; Z0=1.0; Zinf=0.05
phi=0.5; eps_u=5.0; h0=10.0
def s_p(h): return h/(h+h0)
def Z(h): return (1.0+(Zinf/Z0)*h)/(1.0+h)
def k_tot(h): return k0*phi*(1.0+eps_u*h/(h+1.0))+k0*(1.0-phi)/Z(h)
def a_eff(h): return alpha_a*Z(h)/k_tot(h)
def P_eff(h): return P_max*s_p(h)
def Lam_h(h): return a_eff(h)*P_eff(h)/u_star2

def G_h(S,h,n):
    d=d_of_S(S); u=a_eff(h)*P_eff(h)*d/(d+d0)
    return D_of_u(u,n)
def has_psc(h,n):
    Sg=[i/2000.0 for i in range(2001)]
    H=[G_h(s,h,n)-s for s in Sg]; roots=[]
    for i in range(len(Sg)-1):
        if H[i]*H[i+1]<0:
            lo,hi=Sg[i],Sg[i+1]; Hlo=H[i]
            for _ in range(50):
                mid=0.5*(lo+hi); Hm=G_h(mid,h,n)-mid
                if Hlo*Hm<=0: hi=mid
                else: lo=mid; Hlo=Hm
            roots.append(0.5*(lo+hi))
    for r in roots:
        if r>1e-9:
            eps=1e-5; g=(G_h(r+eps,h,n)-G_h(r-eps,h,n))/(2*eps)
            if abs(g)<1.0: return True
    return False

print("%-6s %-12s %-12s %-10s %-12s"%("n","h_low(start)","h_crit(end)","window width","Λ_max(peak)"))
hs=[0.01,0.02,0.03,0.05,0.08,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.8,1.0,
     2.0,3.0,5.0,8.0,10.0,12.0,15.0,20.0,30.0,50.0,80.0,120.0,200.0,300.0]
for n in [0.5,1.0,2.0,3.0]:
    hl=None; hc=None; prev=False
    for h in hs:
        p=has_psc(h,n)
        if not prev and p and hl is None: hl=h
        if prev and not p and hc is None: hc=h
        prev=p
    Lam_peak=max(Lam_h(h) for h in hs)
    w=(hc-hl) if (hl is not None and hc is not None) else None
    print("%-6.1f %-12s %-12s %-10s %-12.2f"%(
        n, "%.3f"%hl if hl else "none", "%.1f"%hc if hc else "none (never quenches)",
        "%.1f"%w if w else "—", Lam_peak))
