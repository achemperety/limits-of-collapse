"""Adversarial test: directed m-cycle, visited-set polynomial G(v) = sum_X gamma_X v^X with
gamma_X = z^|X| mu(P_{m-|X|})/mu(P_m) on cyclic intervals X (path sets), free elsewhere.
Maximize over the free coefficients the minimum, over sample points, of the Rayleigh differences
Delta_ij(v) = d_iG d_jG - G d_ij G (all pairs).  Brandén: real stable => all Delta_ij >= 0 on R^m.
The theorem predicts max-min < 0 for large z."""
import itertools, sys, numpy as np
from scipy.optimize import minimize
m = int(sys.argv[1]); z = float(sys.argv[2]); seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
def mu_path(k, z):
    a, b = 1.0, z
    if k == 0: return 1.0
    for _ in range(k - 1): a, b = b, z * b - a
    return b
subsets = [frozenset(c) for k in range(m + 1) for c in itertools.combinations(range(m), k)]
def is_interval(X):
    if len(X) in (0, m): return True
    return any(all(((s + i) % m) in X for i in range(len(X))) for s in range(m))
known = {X: z ** len(X) * mu_path(m - len(X), z) / mu_path(m, z) for X in subsets if is_interval(X)}
free = [X for X in subsets if not is_interval(X)]
print('m=%d z=%g free coefficients: %d' % (m, z, len(free)))
idx = {X: i for i, X in enumerate(subsets)}
masks = np.array([[1 if t in X else 0 for t in range(m)] for X in subsets])
rng = np.random.default_rng(seed)
# sample points: random, plus points near -1 which carry the certificate
P = np.vstack([rng.uniform(-3, 2, size=(3000, m)), -1 + rng.normal(0, 0.6, size=(3000, m))])
def coeffs(theta):
    c = np.empty(len(subsets))
    for X, val in known.items(): c[idx[X]] = val
    for j, X in enumerate(free): c[idx[X]] = 1.0 + theta[j]
    return c
pairs = list(itertools.combinations(range(m), 2))
def evalG(c, V):  # G and partial derivatives needed, vectorised over rows of V
    mon = np.ones((V.shape[0], len(subsets)))
    for t in range(m):
        mon *= np.where(masks[:, t][None, :] == 1, V[:, t][:, None], 1.0)
    return mon
def min_rayleigh(theta, V=P):
    c = coeffs(theta)
    worst = np.inf
    for (i, j) in pairs:
        # G = G0 + v_i Gi + v_j Gj + v_i v_j Gij, each a polynomial in the other variables
        Vz = V.copy(); Vz[:, i] = 0; Vz[:, j] = 0
        mon = evalG(c, Vz)
        def part(a, b):
            sel = (masks[:, i] == a) & (masks[:, j] == b)
            # monomials in other vars: divide out nothing since v_i=v_j=0 kills; recompute with v_i=v_j=1
            Vo = V.copy(); Vo[:, i] = 1; Vo[:, j] = 1
            mo = evalG(c, Vo)
            return (mo[:, sel] * c[sel][None, :]).sum(1)
        G0, Gi, Gj, Gij = part(0, 0), part(1, 0), part(0, 1), part(1, 1)
        D = Gi * Gj - G0 * Gij
        scale = 1 + np.abs(Gi * Gj) + np.abs(G0 * Gij)
        worst = min(worst, (D / scale).min())
    return worst
best = None
for trial in range(6):
    th0 = rng.normal(0, 0.05, size=len(free)) if trial else np.zeros(len(free))
    res = minimize(lambda th: -min_rayleigh(th), th0, method='Nelder-Mead',
                   options={'maxiter': 4000, 'xatol': 1e-10, 'fatol': 1e-14})
    val = -res.fun
    if best is None or val > best[0]: best = (val, res.x)
    print('trial %d: best min normalized Rayleigh difference = %.3e' % (trial, val), flush=True)
print('RESULT m=%d z=%g  max over free coefficients of min Delta = %.3e  (negative => no stable completion found)' % (m, z, best[0]))
print('free coefficients at optimum (gamma-1):', np.round(best[1], 5))
