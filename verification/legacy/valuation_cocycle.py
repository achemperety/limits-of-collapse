"""Valuation potentials and set potentials (Theorems 10.6 and 10.13).

check(n, arcs) answers two questions for the digraph D = ([n], arcs):
  * valuation potential: is there an integer Lambda(W) with Lambda(W - v) - Lambda(W) = +1 when (W, v) is an
    N-position and -1 when it is a P-position, at every reachable position?  (outcome cocycle)
  * set potential at z0: is there M(W) with M(W - v) = R(W, v) M(W) at every reachable position, with z = z0?
The second test is used only as a filter: a set potential over Q(z) survives evaluation at z0 = 37/11, which is
never a pole, so a failure at z0 is a proof of non-existence. Positive answers are re-checked exactly over Q(z)
by setpot_symbolic.set_potential_exact.
"""
import pickle, sys
from functools import lru_cache
from fractions import Fraction
import dvg


def check(n, arcs, z0=Fraction(37, 11)):
    out = dvg.out_nbrs(n, arcs)

    @lru_cache(maxsize=None)
    def win(W, v):
        return any(not win(W - {v}, t) for t in out[v] if t in W and t != v)

    @lru_cache(maxsize=None)
    def R(W, v):
        s = sum((R(W - {v}, t) for t in out[v] if t in W and t != v), Fraction(0))
        return 1 / (z0 - s)

    P = dvg.reachable_positions(n, arcs)
    # propagate Lambda and M along the edges W -> W - v; every reachable W is reached from V by a play,
    # so processing positions by decreasing |W| visits every edge once
    lam = {frozenset(range(n)): 0}
    M = {frozenset(range(n)): Fraction(1)}
    ok_l = ok_m = True
    for (W, v) in sorted(P, key=lambda p: -len(p[0])):
        if W not in lam:
            continue
        d = 1 if win(W, v) else -1
        W2 = W - {v}
        if W2 in lam:
            ok_l &= (lam[W2] == lam[W] + d)
        else:
            lam[W2] = lam[W] + d
        m2 = M[W] * R(W, v)
        if W2 in M:
            ok_m &= (M[W2] == m2)
        else:
            M[W2] = m2
    return ok_l, ok_m


if __name__ == '__main__':
    res = {}
    for n in (3, 4):
        for arcs in dvg.all_digraphs(n):
            l, m = check(n, arcs)
            key = (n, dvg.scc_symmetric(n, arcs), l, m)
            res[key] = res.get(key, 0) + 1
    for k in sorted(res):
        print('n=%d SCC-sym=%s valuation-potential=%s set-potential(z0)=%s : %d' % (k[0], k[1], k[2], k[3], res[k]))
    print('--- non-SCC-symmetric digraphs with a set potential (n<=4):')
    for n in (3, 4):
        for arcs in dvg.all_digraphs(n):
            if dvg.scc_symmetric(n, arcs):
                continue
            l, m = check(n, arcs)
            if m:
                print(n, arcs)
