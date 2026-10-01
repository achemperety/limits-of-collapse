"""Symbolic checks of the cycle and book theorems (self-energies as indeterminates), and two more predictions."""
from functools import lru_cache
from sage.all__sagemath_singular import PolynomialRing, QQ, FractionField
from setpot_symbolic import set_potential_exact

def gammas(k, arcs, eta, FA, zz):
    out = {v: [w for (a, w) in arcs if a == v] for v in range(k)}
    @lru_cache(maxsize=None)
    def R(W, v):
        s = eta[v]
        for w in out[v]:
            if w in W and w != v: s += R(W - {v}, w)
        return 1 / (zz - s)
    groups = {}
    def ext(path, W, g):
        v = path[-1]; g2 = g * R(W, v)
        groups.setdefault(frozenset(path), set()).add(g2)
        for w in out[v]:
            if w in W and w != v and w not in path: ext(path + [w], W - {v}, g2)
    for s in range(k): ext([s], frozenset(range(k)), FA(1))
    return groups

P = PolynomialRing(QQ, ['z', 'a', 'b', 'e']); FA = FractionField(P)
z, a, b, e = [FA(x) for x in P.gens()]
print('cycles with 2-periodic self-energies (a, b, a, b, ...):')
for m in range(3, 9):
    arcs = [(i, (i + 1) % m) for i in range(m)]
    eta = [a if i % 2 == 0 else b for i in range(m)] if m % 2 == 0 else [a] * m
    g = gammas(m, arcs, eta, FA, z)
    print('   m = %d: set potential for all a%s: %s' % (m, ', b' if m % 2 == 0 else '', all(len(v) == 1 for v in g.values())))
print('books with r pages: eta_y = eta_p = e, eta_x = e + (r-1)/(z-e):')
for r in range(2, 7):
    k = r + 2
    arcs = [(0, 1)] + [(1, p) for p in range(2, k)] + [(p, 0) for p in range(2, k)]
    eta = [e + (r - 1) / (z - e)] + [e] * (r + 1)
    g = gammas(k, arcs, eta, FA, z)
    bad = [X for X, v in g.items() if len(v) > 1]
    eta2 = [e + r / (z - e)] + [e] * (r + 1)
    g2 = gammas(k, arcs, eta2, FA, z)
    print('   r = %d: set potential: %s | with (r)/(z-e) instead: %s' % (r, not bad, all(len(v) == 1 for v in g2.values())))
def sym(u, v): return [(u, v), (v, u)]
blow = [(0, 2), (2, 1), (1, 3), (1, 4), (3, 0), (4, 0)]
print('C4 blow-up (1,1,1,2) with sinks at the vertices of S0 and S1 (predicted True):',
      set_potential_exact(7, blow + [(0, 5), (2, 6)]))
print('C4 blow-up (1,1,1,2) with sinks at S0 only (predicted False):', set_potential_exact(6, blow + [(0, 5)]))
book2 = [(0, 1), (1, 2), (1, 3), (2, 0), (3, 0)]
print('2-page book, sinks at y, p, q; at x a sink and a pendant that has its own sink (predicted True):',
      set_potential_exact(10, book2 + [(1, 4), (2, 5), (3, 6), (0, 7)] + sym(0, 8) + [(8, 9)]))
