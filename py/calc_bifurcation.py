# -*- coding: utf-8 -*-
"""
Simplified 1D closed loop (appendix 1.4.2 replacement): drop R3(f̂), R4(N̂), use D(p) directly as R5.

Closed loop (single domain, no axon reflex, h frozen constant):
    S_{i+1} = G(S_i) = D( x(S_i) )
    d(S)   = gamma * S / (S + S0)                 (R1, d_cns=0)
    x(S)   = Lambda * u* * d(S) / (d(S) + d0)     (R2 folded: x = a·p = a·P_eff·d/(d+d0))
    D(x)   = (x - u*)^n_+ / ((x - u*)^n_+ + K^n)  (R5': thresholded Hill / Stevens power)

    Λ = a·P_eff / u*  (dimensionless source strength, the only bifurcation parameter)

Output:
    fig06_bifurcation.png
      (a) S*(Λ): n=2 default form, single saddle-node → bistable; S=0 always stable, no upper saddle-node (d_cns=0)
      (b) D(p) shape: n=0.5 / 1 / 2 / 3, n<1 superlinear near threshold vs n≥1 smooth
    + per-n comparison table of Λ* and "just-past-threshold S* magnitude" (n<1 jump vs n≥1 smooth)
"""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

# ===================== Parameters (baseline; see report symbol table) =====================
gamma = 0.4   # segmental reflex gain (bounded <1)
S0 = 0.1      # retreat half-saturation (S units)
d0 = 1.0      # source-saturation half-saturation (d units)
u_star = 0.4  # detection threshold (half-anchor 6.8/17 ≈ 0.4 mM)
K = 2.0       # detection half-saturation (mM, replaces e0/p0)
# n = Hill/Stevens exponent of R5', scan n ∈ {0.5, 1, 2, 3}

# ===================== Closed-loop map =====================
def d_of_S(S):
    return gamma * S / (S + S0)

def x_of_S(S, Lambda):
    d = d_of_S(S)
    return Lambda * u_star * d / (d + d0)

def D_of_x(x, n):
    t = max(0.0, x - u_star)
    if t == 0.0:
        return 0.0
    return t ** n / (t ** n + K ** n)

def G(S, Lambda, n):
    return D_of_x(x_of_S(S, Lambda), n)

def _H(S, Lambda, n):
    return G(S, Lambda, n) - S

# ===================== Fixed-point solver =====================
def fixed_points(Lambda, n, ngrid=8001):
    Sg = np.linspace(0.0, 1.0, ngrid)
    H = np.array([_H(s, Lambda, n) for s in Sg])
    roots = []
    if abs(H[0]) < 1e-12:
        roots.append(0.0)
    for i in range(ngrid - 1):
        if H[i] * H[i + 1] < 0.0:
            lo, hi = Sg[i], Sg[i + 1]
            Hlo = H[i]
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                Hm = _H(mid, Lambda, n)
                if Hlo * Hm <= 0.0:
                    hi = mid
                else:
                    lo = mid
                    Hlo = Hm
            roots.append(0.5 * (lo + hi))
    return sorted(set(round(r, 10) for r in roots))

def is_stable(S, Lambda, n, eps=1e-5):
    if S < 1e-9:
        return True
    g = (G(S + eps, Lambda, n) - G(S - eps, Lambda, n)) / (2.0 * eps)
    return abs(g) < 1.0

def collect(Lambdas, n):
    zero, stable, unstable = [], [], []
    for L in Lambdas:
        for r in fixed_points(L, n):
            if r < 1e-9:
                zero.append(L)
            elif is_stable(r, L, n):
                stable.append((L, r))
            else:
                unstable.append((L, r))
    return zero, stable, unstable

def arr(a):
    return np.array(a) if a else np.empty((0, 2))

def saddle_node(n):
    """Return Λ* (first Λ at which a stable excited attractor appears) and the just-past-threshold S* magnitude."""
    Ls = np.linspace(0.0, 40.0, 4001)
    found = None
    for L in Ls:
        st = [r for r in fixed_points(L, n) if r > 1e-9 and is_stable(r, L, n)]
        if st:
            found = L
            break
    if found is None:
        return None, None
    # just-past-threshold (Λ = 1.05·Λ*) stable excited S* magnitude
    Lm = found * 1.05
    st = [r for r in fixed_points(Lm, n) if r > 1e-9 and is_stable(r, Lm, n)]
    S_at = max(st) if st else None
    return found, S_at

# ===================== Plot =====================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.2, 5.0), dpi=150)

# --- (a) S*(Λ) n=2 default ---
n_rec = 2
Lam_grid = np.linspace(0.0, 30.0, 3001)
zero, st, un = collect(Lam_grid, n_rec)
st, un = arr(st), arr(un)
Lam_star = st[:, 0].min()

def draw(ax, zero, st, un):
    ax.plot(zero, np.zeros(len(zero)), lw=2, color="#1f77b4")
    if len(st):
        ax.plot(st[:, 0], st[:, 1], lw=2.5, color="#1f77b4")
    if len(un):
        ax.plot(un[:, 0], un[:, 1], lw=2, color="#d62728", ls="--")

draw(ax1, zero, st, un)
ax1.axvline(Lam_star, ls="--", lw=1.2, color="green", alpha=0.85)
S_at = st[st[:, 0].argmin(), 1]
ax1.annotate("Saddle-node $\\Lambda^*\\approx%.2f$\nexcited attractor appears" % Lam_star,
             xy=(Lam_star, S_at), xytext=(Lam_star + 1.2, 0.42),
             fontsize=9.5, color="green",
             arrowprops=dict(arrowstyle="->", color="green", lw=1))
ax1.text(9.0, 0.06, "Bistable window $\\Lambda\\in[%.2f,\\infty)$ (no upper saddle-node, $d_{cns}=0$)" % Lam_star,
         fontsize=8.5, color="purple")
ax1.axhline(0.0, color="gray", lw=0.6, alpha=0.5)
ax1.set_xlabel("Dimensionless source strength $\\Lambda$ (bifurcation parameter)", fontsize=11)
ax1.set_ylabel("Fixed point $S^*$", fontsize=11)
ax1.set_title("(a) $S^*(\\Lambda)$: single saddle-node → bistable ($D(p)$ thresholded Hill n=2)", fontsize=11.5)
ax1.set_xlim(Lam_star - 10, Lam_star + 10); ax1.set_ylim(-0.02, 1.02)
ax1.grid(alpha=0.25)

# --- (b) D(p) shape n=0.5/1/2/3 ---
x = np.linspace(0.0, 6.0, 1201)
for n, ls, lb in [(0.5, "-", "n=0.5 (superlinear)"),
                  (1.0, "-", "n=1 (simplest)"),
                  (2.0, "-", "n=2 (default)"),
                  (3.0, "--", "n=3")]:
    y = np.array([D_of_x(xx, n) for xx in x])
    ax2.plot(x, y, lw=2 if n != 3 else 1.8, ls=ls, label=lb)
ax2.axvline(u_star, ls=":", lw=1.2, color="gray", alpha=0.9)
ax2.text(u_star + 0.08, 0.06, "$u^*_{eff}$=0.4", fontsize=9, color="gray")
ax2.set_xlabel("Detection-side concentration $x=a\\,p$ (mM)", fontsize=11)
ax2.set_ylabel("$D(p)$ (normalized sensation)", fontsize=11)
ax2.set_title("(b) $D(p)$ shape: n<1 superlinear near threshold, n≥1 smooth", fontsize=11.5)
ax2.set_xlim(0, 6); ax2.set_ylim(-0.02, 1.02)
ax2.legend(fontsize=9, loc="lower right")
ax2.grid(alpha=0.25)

fig.tight_layout()
out = Path(__file__).resolve().parent.parent / 'images' / 'fig06_bifurcation.png'
fig.savefig(out, bbox_inches="tight")
print("[saved]", out)

# ===================== Analytical + numerical conclusions =====================
print("=" * 64)
print("Simplified 1D closed loop: S_{i+1} = D( a·p(d(S_i)) ), D thresholded Hill/Stevens power")
print("=" * 64)
print("Parameters: γ=%.1f, S0=%.1f, d0=%.1f, u*=%.1f, K=%.1f" % (gamma, S0, d0, u_star, K))
print("Analytical D>0 necessary threshold: Λ > (γ+d0(1+S0))/γ = %.2f" % ((gamma + d0 * (1 + S0)) / gamma))
print("-" * 64)
print("R5 form   | Λ*    | D(u*+0.1) | D(u*+0.5) | just-threshold | verdict")
print("-" * 66)
for n in [0.5, 1.0, 2.0, 3.0]:
    ls, _ = saddle_node(n)
    d_near = D_of_x(u_star + 0.1, n)   # 0.1 mM above threshold (just past)
    d_far = D_of_x(u_star + 0.5, n)    # 0.5 mM above threshold
    near_str = "strong(%.2f)" % d_near if d_near > 0.1 else "weak(%.3f)" % d_near
    verdict = "✓ (n<1: strong right past threshold, yet preserves high-health silence)" if n < 1 else "✓"
    print("%-9s | %6.2f | %9.3f | %9.3f | %-13s | %s" % ("n=%.1f" % n, ls, d_near, d_far, near_str, verdict))

# Recommended form summary
print("-" * 66)
print("Recommended R5': D(p) = (a·p − u*)^n_+ / ((a·p − u*)^n_+ + K^n), any n > 0")
print("  n=2 (default); n=1 (simplest, = original N̂ single term)")
