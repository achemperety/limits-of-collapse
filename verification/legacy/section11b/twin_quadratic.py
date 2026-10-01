"""Twins in a single-class blow-up: the rotation through t_j at t_j and at x_{m-1} give an equation for
v_j = z - eta(t_j) whose coefficients depend on the twins only through symmetric sums. Check it and its roots."""
from functools import lru_cache
from sage.all__sagemath_singular import PolynomialRing, QQ, FractionField
for m, r in [(3, 2), (3, 3), (4, 2), (5, 2)]:
    P = PolynomialRing(QQ, ['z'] + ['v%d' % j for j in range(r)] + ['x%d' % i for i in range(1, m)])
    F = FractionField(P); z = F(P.gen(0)); vs = [F(P.gen(1 + j)) for j in range(r)]
    xs = [None] + [F(P.gen(1 + r + i - 1)) for i in range(1, m)]     # these are weights w_i = z - eta
    tw = list(range(r)); X = {i: r + i - 1 for i in range(1, m)}; n = r + m - 1
    arcs = [(t, X[1]) for t in tw] + [(X[m - 1], t) for t in tw] + [(X[i], X[i + 1]) for i in range(1, m - 1)]
    wt = {t: vs[t] for t in tw}; wt.update({X[i]: xs[i] for i in range(1, m)})
    out = {v: [b for a, b in arcs if a == v] for v in range(n)}
    @lru_cache(maxsize=None)
    def R(W, v):
        s = F(0)
        for u in out[v]:
            if u in W and u != v: s += R(W - {v}, u)
        return 1 / (wt[v] - s)
    def G(path):
        W = frozenset(range(n)); g = F(1)
        for v in path: g *= R(W, v); W = W - {v}
        return g
    j = 0
    g1 = G([tw[j]] + [X[i] for i in range(1, m)])            # rotation at t_j
    g2 = G([X[m - 1], tw[j]] + [X[i] for i in range(1, m - 1)])  # rotation at x_{m-1} through t_j
    num = (1 / g1 - 1 / g2).numerator()
    Pv = PolynomialRing(FractionField(PolynomialRing(QQ, [str(g) for g in P.gens() if str(g) != 'v0'])), 'v0')
    print('m=%d r=%d: degree of the equation in v_j:' % (m, r), Pv(num).degree())
