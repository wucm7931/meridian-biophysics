# -*- coding: utf-8 -*-
"""Derive the asymptotic form of f(Λ*-Λ): saddle-node ghost / critical slowing N≈C/√δ, verifying C=π/√(H_Λ|H_SS|/2)."""
import sys, math
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

u_star=0.4; gamma=0.4; S0=0.1; d0=1.0; e0=2.0; n=2.0; S_seed=0.8

def d_of_S(S): return gamma*S/(S+S0)
def m_of_S(S): return d_of_S(S)/(d_of_S(S)+d0)
def u_of_S(S, Lam): return Lam*u_star * m_of_S(S)
def D_of_u(u):
    x=max(0.0,u-u_star); return x**n/(x**n+e0**n)
def G(S, Lam): return D_of_u(u_of_S(S, Lam))

# ---- 1) saddle-node Λ* ----
def has_stable_nonzero(Lam):
    Sg=[i/4000.0 for i in range(4001)]
    Hl=[G(s,Lam)-s for s in Sg]; roots=[]
    for i in range(len(Sg)-1):
        if Hl[i]*Hl[i+1]<0:
            lo,hi=Sg[i],Sg[i+1]; Hlo=Hl[i]
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

# ---- 2) find S* at Λ* (G(S*)=S*, tangency point G'(S*)=1) ----
Sg=[i/200000.0 for i in range(200001)]
Hs=[G(s,Lam_star)-s for s in Sg]
S_star=None
for i in range(len(Sg)-1):
    if Hs[i]*Hs[i+1]<0:
        l,h2=Sg[i],Sg[i+1]; Hl=Hs[i]
        for _ in range(70):
            m=0.5*(l+h2); Hm=G(m,Lam_star)-m
            if Hl*Hm<=0: h2=m
            else: l=m; Hl=Hm
        S_star=0.5*(l+h2); break

# ---- 3) H_Λ and H_SS (finite differences) ----
epsL=1e-5; epsS=1e-5
H_Lam = (G(S_star, Lam_star+epsL) - G(S_star, Lam_star-epsL))/(2*epsL)
H_SS  = (G(S_star+epsS, Lam_star) - 2*G(S_star, Lam_star) + G(S_star-epsS, Lam_star))/(epsS**2)

C = math.pi / math.sqrt(H_Lam * abs(H_SS) / 2.0)

print("Λ* = %.5f"%Lam_star)
print("S* = %.5f"%S_star)
print("G'(S*) = %.6f  (should be 1)"%((G(S_star+1e-5,Lam_star)-G(S_star-1e-5,Lam_star))/(2e-5)))
print("H_Λ = dG/dΛ = %.6f"%H_Lam)
print("H_SS = d2G/dS2 = %.6f"%H_SS)
print("C = pi/sqrt(H_Λ|H_SS|/2) = %.4f"%C)
print()

# ---- 4) numerical segment count vs prediction C/√δ ----
def decay_segments(Lam):
    S=S_seed
    for k in range(1000000):
        if u_of_S(S,Lam) < u_star: return k
        S=G(S,Lam)
    return 1000000

print("δ=Λ*-Λ    numerical segments N     C/√δ       N·√δ (should →C)")
print("-"*58)
for dlt in [0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]:
    N=decay_segments(Lam_star-dlt)
    pred=C/math.sqrt(dlt)
    print("%-10.4f %-11d %-10.2f %-10.3f"%(dlt, N, pred, N*math.sqrt(dlt)))
