import itertools, numpy as np, sys
sys.argv = ['x', '5', '8', '1', '800', '1']
exec(open('adversarial_ctrl.py').read().split('best = None')[0])
# exact completion for the bidirected 5-cycle: gamma_X = z^|X| mu(C_m - X)/mu(C_m)
def mu_graph_cycle_minus(X):
    # C_m minus X is a disjoint union of paths; multiply mu of the runs of C\X
    rest = [t for t in range(m) if t not in X]
    if len(rest) == m: return mu_cyc(m, z)
    # runs of consecutive vertices in cyclic order
    runs, cur = [], 0
    start = [t for t in rest if (t - 1) % m in X][0]
    t = start; L = 0; out = 1.0
    for k in range(m):
        u = (start + k) % m
        if u in X:
            if L: out *= mu_path(L, z); L = 0
        else: L += 1
    if L: out *= mu_path(L, z)
    return out
theta = np.array([z ** len(X) * mu_graph_cycle_minus(X) / mu_cyc(m, z) - 1.0 for X in free])
print('exact bidirected-cycle completion: min normalized Delta =', min_rayleigh(theta))
print('theta:', np.round(theta, 4))
