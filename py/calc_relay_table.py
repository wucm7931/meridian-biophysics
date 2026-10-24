# -*- coding: utf-8 -*-
"""Point ⑤ five-domain relay table (n=2, γ=0.4 unified convention): self-sustaining branch, S converges from seed 0.8 toward S*."""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

# anchored parameters (n=2, γ=0.4 unified)
gamma=0.4; S0=0.1; d0=1.0; u_star=0.4; K=2.0; n=2.0
P_eff=16.0; a=0.40; S_seed=0.8

def d_of_S(S): return gamma*S/(S+S0)
def p_of_d(d): return P_eff*d/(d+d0)
def D_of_u(u):
    x=max(0.0, u-u_star); return x**n/(x**n+K**n)

Lam = a*P_eff/u_star
print("P_eff=%.1f, a=%.2f, u*_eff=%.1f -> Λ=%.3f  (Λ*=15.892)"%(P_eff,a,u_star,Lam))
print()

# find stable nonzero fixed point S* (G(S*)=S*, |G'|<1)
def G(S):
    d=d_of_S(S); p=p_of_d(d); u=a*p; return D_of_u(u)

Sg=[i/100000.0 for i in range(100001)]
roots=[]
Hs=[G(s)-s for s in Sg]
for i in range(len(Sg)-1):
    if Hs[i]*Hs[i+1]<0:
        lo,hi=Sg[i],Sg[i+1]; Hlo=Hs[i]
        for _ in range(60):
            mid=0.5*(lo+hi); Hm=G(mid)-mid
            if Hlo*Hm<=0: hi=mid
            else: lo=mid; Hlo=Hm
        roots.append(0.5*(lo+hi))
S_star=max(r for r in roots if r>1e-9)
print("stable nonzero fixed point S* = %.4f"%S_star)
print()

# five-domain relay
print("%-10s %-14s %-12s %-12s %-12s"%("domain","u(detection,mM)","p(source peak,mM)","S(sensation)","d(retreat)"))
S=S_seed
for j,name in enumerate(["i(seed)","i+1","i+2","i+3","i+4"]):
    d=d_of_S(S); p=p_of_d(d); u=a*p; Snext=D_of_u(u)
    print("%-10s %-14.3f %-12.3f %-12.3f %-12.3f"%(name,u,p,S,d))
    S=Snext
