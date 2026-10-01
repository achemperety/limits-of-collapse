"""Adversarial test (fast version).  Directed m-cycle; visited-set polynomial G(v)=sum_X gamma_X v^X with
gamma_X known on cyclic intervals, free elsewhere. Maximize over free coefficients the minimum over sample
points of the normalized Rayleigh differences. Theorem predicts a negative optimum for large z."""
import itertools, sys, numpy as np
from scipy.optimize import minimize
m = int(sys.argv[1]); z = float(sys.argv[2]); seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
npts = int(sys.argv[4]) if len(sys.argv) > 4 else 1500
def mu_path(k, z):
    a, b = 1.0, z
    if k == 0: return 1.0
    for _ in range(k - 1): a, b = b, z * b - a
    return b
subsets = [frozenset(c) for k in range(m + 1) for c in itertools.combinations(range(m), k)]
def is_interval(X):
    if len(X) in (0, m): return True
    return any(all(((s + i) % m) in X for i in range(len(X))) for s in range(m))
def mu_cyc(m, z): return mu_path(m, z) - mu_path(m - 2, z)
SYM = int(sys.argv[5]) if len(sys.argv) > 5 else 0
den = mu_cyc(m, z) if SYM else mu_path(m, z)
known = {X: (1.0 if len(X) == 0 else z ** len(X) * (mu_path(m - len(X), z) if len(X) < m else 1.0) / den) for X in subsets if is_interval(X)}
free = [X for X in subsets if not is_interval(X)]
N = len(subsets); idx = {X: i for i, X in enumerate(subsets)}
rng = np.random.default_rng(seed)
P = np.vstack([rng.uniform(-3, 2, size=(npts, m)), -1 + rng.normal(0, 0.6, size=(npts, m))])
pairs = list(itertools.combinations(range(m), 2))
blocks = []
for (i, j) in pairs:
    mats = []
    for a, b in ((0, 0), (1, 0), (0, 1), (1, 1)):
        Mab = np.zeros((P.shape[0], N))
        for X in subsets:
            if ((i in X) == bool(a)) and ((j in X) == bool(b)):
                rest = [t for t in X if t not in (i, j)]
                Mab[:, idx[X]] = np.prod(P[:, rest], axis=1) if rest else 1.0
        mats.append(Mab)
    blocks.append(mats)
base = np.zeros(N)
for X, val in known.items(): base[idx[X]] = val
Fcols = [idx[X] for X in free]
def coeffs(theta):
    c = base.copy(); c[Fcols] = 1.0 + theta; return c
def min_rayleigh(theta):
    theta = np.asarray(theta, dtype=float)
    c = coeffs(theta); worst = np.inf
    for M0, Mi, Mj, Mij in blocks:
        g0, gi, gj, gij = M0 @ c, Mi @ c, Mj @ c, Mij @ c
        D = gi * gj - g0 * gij
        worst = min(worst, (D / (1 + np.abs(gi * gj) + np.abs(g0 * gij))).min())
    return worst
best = None
if not free:
    print('m=%d z=%g free=0  min normalized Delta = %.3e' % (m, z, min_rayleigh(np.zeros(0)))); sys.exit()
for trial in range(8):
    th0 = rng.normal(0, 0.2 if trial > 3 else 0.02, size=len(free)) if trial else np.zeros(len(free))
    res = minimize(lambda th: -min_rayleigh(th), th0, method='Nelder-Mead',
                   options={'maxiter': 20000, 'maxfev': 20000, 'xatol': 1e-12, 'fatol': 1e-15})
    val = -res.fun
    if best is None or val > best[0]: best = (val, res.x)
print('m=%d z=%g free=%d  max_free min_points normalized Delta = %.3e   (z^-2 = %.1e)' % (m, z, len(free), best[0], z ** -2))
