"""Theorem 11.11: the explicit Gamma of every path of the 3-page book with independent self-energies."""
from functools import lru_cache
from sage.all__sagemath_singular import PolynomialRing, QQ, FractionField
P = PolynomialRing(QQ, ['z', 'ex', 'ey', 'e2', 'e3', 'e4']); FA = FractionField(P)
z, ex, ey, e2, e3, e4 = [FA(g) for g in P.gens()]
arcs = [(0, 1), (1, 2), (1, 3), (1, 4), (2, 0), (3, 0), (4, 0)]
eta = {0: ex, 1: ey, 2: e2, 3: e3, 4: e4}
out = {v: [w for a, w in arcs if a == v] for v in range(5)}
@lru_cache(maxsize=None)
def R(W, v):
    s = eta[v]
    for w in out[v]:
        if w in W and w != v: s += R(W - {v}, w)
    return 1 / (z - s)
def G(path):
    W = frozenset(range(5)); g = FA(1)
    for v in path: g *= R(W, v); W = W - {v}
    return g
zx, zy = z - ex, z - ey
zp = {2: z - e2, 3: z - e3, 4: z - e4}
S = sum(1 / zp[j] for j in zp); Tt = sum(zx / (zp[j] * zx - 1) for j in zp)
ok = True
for i in zp:
    Si = S - 1 / zp[i]
    Di = zp[i] * (zx * (zy - Si) - 1) - (zy - Si)
    ok &= G((0, 1, i)) == 1 / ((zx * (zy - S) - 1) * zp[i])
    ok &= G((1, i, 0)) == 1 / ((zy - Tt) * (zp[i] * zx - 1))
    ok &= G((i, 0, 1)) == 1 / Di
    for j in zp:
        if j != i: ok &= G((i, 0, 1, j)) == 1 / (Di * zp[j])
print('explicit Gamma formulas of Theorem 11.11 hold for the 3-page book:', ok)
