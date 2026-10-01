"""Set potentials decided exactly over Q(z) (Theorem 10.8 pruning and the census of Theorem 10.13).

set_potential_exact(n, arcs) propagates M from M(V) = 1 along the reachable positions, M(W - v) = R(W, v) M(W),
and reports whether every position is consistent. The ratios R(W, v) lie in Q(z), so a nowhere-zero solution in
any extension field of Q(z) exists exactly when this test succeeds.

Run as a script: the census of set potentials on non-SCC-symmetric digraphs with 3, 4 and 5 vertices
(expected: 86 flagged at z0 = 37/11, all 86 confirmed exactly; 1 + 6 + 79 by size).
"""
import pickle
import dvg
from valuation_cocycle import check


def set_potential_exact(n, arcs):
    R = dvg.make_resolvent(n, arcs)          # exact over Q(z)
    P = dvg.reachable_positions(n, arcs)
    M = {frozenset(range(n)): dvg.Kz(1)}
    for (W, v) in sorted(P, key=lambda p: -len(p[0])):
        if W not in M:
            continue
        m2 = M[W] * R(W, v)
        W2 = W - {v}
        if W2 in M:
            if M[W2] != m2:
                return False
        else:
            M[W2] = m2
    return True


if __name__ == '__main__':
    flag = []
    for n in (3, 4):
        for arcs in dvg.all_digraphs(n):
            if not dvg.scc_symmetric(n, arcs) and check(n, arcs)[1]:
                flag.append((n, arcs))
    for arcs in pickle.load(open('digraphs5.pkl', 'rb')):
        if not dvg.scc_symmetric(5, arcs) and check(5, arcs)[1]:
            flag.append((5, arcs))
    ok = sum(set_potential_exact(n, a) for n, a in flag)
    by_n = {k: sum(1 for n, a in flag if n == k) for k in (3, 4, 5)}
    print('set potentials flagged at z0:', len(flag), by_n, '| confirmed exactly over Q(z):', ok)
