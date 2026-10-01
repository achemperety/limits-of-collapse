"""For every strongly connected, non-symmetric digraph C on k vertices, find all self-energies eta (indeterminates
over Q(z)) for which the game on C with self-energies eta has a set potential. Output: the prime components of the
solution variety (Groebner basis over Q(z) after saturation by all resolvent denominators)."""
import sys, itertools, pickle
from functools import lru_cache
from sage.all__sagemath_singular import PolynomialRing, QQ, FractionField
import dvg

k = int(sys.argv[1]) if __name__ == '__main__' else 5
Rz = PolynomialRing(QQ, 'z'); Kz = FractionField(Rz); zz = Kz.gen()
A = PolynomialRing(Kz, ['e%d' % i for i in range(k)], order='degrevlex')
E = A.gens()
FA = FractionField(A)


def strongly_connected(n, arcs):
    r = dvg.reach_matrix(n, arcs)
    return all(r[a][b] for a in range(n) for b in range(n))


def analyse(arcs):
    out = {v: [w for (a, w) in arcs if a == v] for v in range(k)}

    @lru_cache(maxsize=None)
    def R(W, v):
        s = FA(E[v])
        for w in out[v]:
            if w in W and w != v:
                s += R(W - {v}, w)
        return 1 / (FA(zz) - s)
    # all directed paths, grouped by vertex set; Gamma of each
    groups = {}
    V = frozenset(range(k))

    def ext(path, W, g):
        v = path[-1]
        g2 = g * R(W, v)
        groups.setdefault(frozenset(path), []).append((tuple(path), g2))
        for w in out[v]:
            if w in W and w != v and w not in path:
                ext(path + [w], W - {v}, g2)
    for s in range(k):
        ext([s], V, FA(1))
    eqs, dens = [], set()
    for X, lst in groups.items():
        g0 = lst[0][1]
        for p, g in lst[1:]:
            d = g - g0
            if d != 0:
                eqs.append(d.numerator())
    # denominators: every resolvent denominator must be nonzero
    for key in list(groups):
        for p, g in groups[key]:
            dens.add(g.denominator())
            dens.add(g.numerator())
    if not eqs:
        return 'free', None
    I = A.ideal(eqs)
    for d in dens:
        if d.is_constant():
            continue
        I = I.saturation(A.ideal([d]))[0]
        if I.is_one():
            return 'none', None
    return 'variety', I


if __name__ == '__main__':
    pairs = [(a, b) for a in range(k) for b in range(k) if a != b]
    res = {}
    seen = set()
    for arcs in dvg.all_digraphs(k):
        if not strongly_connected(k, arcs) or dvg.scc_symmetric(k, arcs):
            continue
        kind, I = analyse(arcs)
        if kind == 'none':
            res.setdefault('none', []).append(arcs)
            continue
        if kind == 'free':
            print('FREE (no condition):', arcs, flush=True)
            continue
        comps = I.minimal_associated_primes()
        print('arcs', arcs, 'dim', I.dimension(), flush=True)
        for P in comps:
            print('    component dim', P.dimension(), ':', [str(g) for g in P.gens()], flush=True)
    print('strongly connected non-symmetric digraphs with no admissible eta:', len(res.get('none', [])))
