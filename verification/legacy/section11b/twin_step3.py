"""Theorem 11.22, Steps 2 and 3, checked with all weights symbolic: the products Gamma of the rotations of Z_j at t_j
and at x_1 are the stated continuants, and they agree exactly when K(L')/v_j = K(L) - K(L', y)."""
from functools import lru_cache
from sage.all__sagemath_singular import PolynomialRing, QQ, FractionField

def cont(seq):
    a, b = 0, 1                      # K_{-1} = 0, K_0 = 1
    for s in seq:
        a, b = b, s * b - a
    return b

for m, r in [(3, 2), (3, 3), (4, 2), (4, 3), (5, 2), (6, 2)]:
    P = PolynomialRing(QQ, ['v%d' % j for j in range(r)] + ['w%d' % i for i in range(1, m)])
    F = FractionField(P); vs = [F(P.gen(j)) for j in range(r)]; ws = [None] + [F(P.gen(r + i - 1)) for i in range(1, m)]
    tw = list(range(r)); X = {i: r + i - 1 for i in range(1, m)}; n = r + m - 1
    arcs = [(t, X[1]) for t in tw] + [(X[m - 1], t) for t in tw] + [(X[i], X[i + 1]) for i in range(1, m - 1)]
    wt = {t: vs[t] for t in tw}; wt.update({X[i]: ws[i] for i in range(1, m)})
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
    Sig = sum(1 / v for v in vs); Sj = Sig - 1 / vs[j]
    L = [ws[i] for i in range(1, m - 1)]; Lp = L[1:]; y = ws[m - 1] - Sig
    ok1 = G([tw[j]] + [X[i] for i in range(1, m)]) == 1 / cont([vs[j]] + L + [ws[m - 1] - Sj])
    ok2 = G([X[i] for i in range(1, m)] + [tw[j]]) == 1 / (vs[j] * cont(L + [y]))
    ok3 = G([tw[j]] + [X[i] for i in range(1, m)] + [tw[1]]) == 1 / (vs[1] * cont([vs[j]] + L + [ws[m - 1] - Sj]))
    lhs = 1 / G([tw[j]] + [X[i] for i in range(1, m)]) - 1 / G([X[i] for i in range(1, m)] + [tw[j]])
    rhs = cont(Lp) / vs[j] - cont(L) + cont(Lp + [y])
    ok4 = (lhs + rhs == 0) or (lhs - rhs == 0) or (lhs / rhs in QQ)
    print('m=%d r=%d  Step 2 formulas: %s %s %s | Step 3 equation equivalent: %s' % (m, r, ok1, ok2, ok3, ok4), flush=True)

# Step 4 with equal twins: rotations at x_a give 1/K_a (coupling r between w_{m-1} and w), the rotation at t_j gives 1/K_t
def contc(seq, coup):
    a, b = 0, 1
    for s, c in zip(seq, coup):
        a, b = b, s * b - c * a
    return b

for m, r in [(3, 2), (4, 2), (4, 3), (5, 2), (5, 3), (6, 2), (7, 2)]:
    P = PolynomialRing(QQ, ['w'] + ['w%d' % i for i in range(1, m)])
    F = FractionField(P); w = F(P.gen(0)); ws = [None] + [F(P.gen(i)) for i in range(1, m)]
    tw = list(range(r)); X = {i: r + i - 1 for i in range(1, m)}; n = r + m - 1
    arcs = [(t, X[1]) for t in tw] + [(X[m - 1], t) for t in tw] + [(X[i], X[i + 1]) for i in range(1, m - 1)]
    wt = {t: w for t in tw}; wt.update({X[i]: ws[i] for i in range(1, m)})
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
    u = [w] + [ws[i] for i in range(1, m)]
    beta = [r] + [1] * (m - 1)                       # beta_i couples u_{i-1} and u_i
    ok = True
    for a in range(1, m):
        seq = [u[(a + k) % m] for k in range(m)]
        coup = [0] + [beta[(a + k) % m] for k in range(1, m)]
        Ka = contc(seq, coup)
        ok &= G([X[i] for i in range(a, m)] + [tw[0]] + [X[i] for i in range(1, a)]) == 1 / Ka
    Kt = contc([w] + [ws[i] for i in range(1, m - 1)] + [ws[m - 1] - (r - 1) / w], [0] + [1] * (m - 1))
    ok &= G([tw[0]] + [X[i] for i in range(1, m)]) == 1 / Kt
    print('m=%d r=%d  Step 4: rotations at x_a are 1/K_a and at t_j is 1/K_t: %s' % (m, r, ok), flush=True)
