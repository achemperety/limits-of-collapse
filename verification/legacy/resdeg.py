"""Reduced degree of the exact DVG root resolvent R(V,s) as a rational function, for bounded-treewidth families."""
import sys, random, itertools
from functools import lru_cache
from sage.all__sagemath_singular import PolynomialRing, QQ, FractionField
Rz = PolynomialRing(QQ, 'z'); z = Rz.gen(); K = FractionField(Rz)
sys.setrecursionlimit(100000)

def resolvent_degree(n, arcs, s):
    out = {v: [] for v in range(n)}
    for a, b in arcs: out[a].append(b)
    @lru_cache(maxsize=None)
    def R(W, v):
        acc = K(0)
        for t in out[v]:
            if t in W and t != v: acc += R(W - {v}, t)
        return 1 / (K(z) - acc)
    r = R(frozenset(range(n)), s)
    npos = R.cache_info().currsize
    return r.denominator().degree(), npos

def ladder(m, orient):
    # 2 x m ladder: rails 0..m-1 (top) and m..2m-1 (bottom); rungs i <-> m+i
    arcs = []
    for i in range(m - 1):
        arcs += [(i, i + 1), (m + i + 1, m + i)] if orient == 'opp' else [(i, i + 1), (m + i, m + i + 1)]
    for i in range(m):
        arcs += [(i, m + i), (m + i, i)]
    return 2 * m, arcs

def cyc_chords(n, rng):
    arcs = [(i, (i + 1) % n) for i in range(n)]
    for _ in range(n // 3):
        a = rng.randrange(n); b = (a + rng.randrange(2, n - 1)) % n
        arcs.append((a, b))
    return n, sorted(set(arcs))

def sp_digraph(n, rng):
    # random 2-tree then random orientation of each edge (one-way or two-way)
    edges = {(0, 1)}
    for v in range(2, n):
        a, b = rng.choice(sorted(edges))
        edges |= {(min(a, v), max(a, v)), (min(b, v), max(b, v))}
    arcs = []
    for a, b in edges:
        r = rng.random()
        if r < 0.4: arcs.append((a, b))
        elif r < 0.8: arcs.append((b, a))
        else: arcs += [(a, b), (b, a)]
    return n, arcs

rng = random.Random(1)
print('ladder (rails opposite, rungs symmetric):')
for m in range(2, 8):
    n, arcs = ladder(m, 'opp'); d, p = resolvent_degree(n, arcs, 0)
    print('  n=%2d  deg R = %4d   positions = %d' % (n, d, p), flush=True)
print('ladder (rails same direction, rungs symmetric):')
for m in range(2, 8):
    n, arcs = ladder(m, 'same'); d, p = resolvent_degree(n, arcs, 0)
    print('  n=%2d  deg R = %4d   positions = %d' % (n, d, p), flush=True)
print('random treewidth-2 digraphs (max over starts):')
for n in range(6, 15, 2):
    best = 0
    for trial in range(3):
        nn, arcs = sp_digraph(n, rng)
        best = max(best, max(resolvent_degree(nn, arcs, s)[0] for s in range(nn)))
    print('  n=%2d  max deg R = %4d' % (n, best), flush=True)
