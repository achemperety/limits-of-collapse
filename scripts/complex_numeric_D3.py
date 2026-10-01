"""Theorem 11.6 for m = 3: the potential identities with the numerical roots of Y_3 (50 digits)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
import time, sys
import mpmath as mp
from dvgcollapse.dvg import Digraph
from dvgcollapse import groebner_saturation as gs
from dvgcollapse.matching import mu_weighted
m = 3
D = Digraph(2 * m + 1, [(i, (i + 1) % m) for i in range(m)])     # C_3 plus 4 isolated vertices
U = [(c, f) for c in range(m) for f in range(m, 2 * m + 1)]
# numerical check at real z = 5 with 50 digits
mp.mp.dps = 50
for zv in (mp.mpf(5), mp.mpf('0.7'), mp.mpc('1.3', '0.4')):
    coeffs = [zv**2 - 2, -zv * (zv**2 - 2), -3 * (zv**2 - 1), -6 * zv, -6]
    roots = mp.polyroots(coeffs, maxsteps=200, extraprec=200)
    a = {0: 0, 1: 0, 2: 0}; a.update({3 + i: roots[i] for i in range(4)})
    V = frozenset(range(7))
    def M(S): return mu_weighted(S, U, a, one=mp.mpf(1))
    def mu_p(k):
        x, y = mp.mpf(1), zv
        if k == 0: return x
        for _ in range(k - 1): x, y = y, zv * y - x
        return y
    worst = 0
    F = frozenset(range(3, 7))
    for s in range(3):
        for L in range(1, 4):
            I = [(s + i) % 3 for i in range(L)]
            r = M(frozenset(I[1:]) | F) / M(frozenset(I) | F) - mu_p(L - 1) / mu_p(L)
            worst = max(worst, abs(r))
    for f in range(3, 7):
        worst = max(worst, abs(M(V - {f}) / M(V) - 1 / zv))
    print('z = %s: env weights %s' % (mp.nstr(zv, 4), [mp.nstr(r, 6) for r in roots]))
    print('     max |identity residual| = %s' % mp.nstr(worst, 3))
