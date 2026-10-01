"""Potentials with a frozen environment: U on V+F, M(W) = mu_a(U[W u F]); constraints only from positions of D."""
import itertools, sys
from sage.all__sagemath_singular import PolynomialRing
import dvg

def has_env_potential(n, arcs, nf, Uedges):
    N = n + nf
    A = PolynomialRing(dvg.Kz, ['a%d' % i for i in range(N)], order='degrevlex')
    al = A.gens(); memo = {}
    F = frozenset(range(n, N))
    P = dvg.reachable_positions(n, arcs)
    R = dvg.make_resolvent(n, arcs)
    eqs = []; Ws = set()
    for (W, v) in P:
        eqs.append(dvg.mu_weighted((W - {v}) | F, Uedges, al, A, memo) - A(R(W, v)) * dvg.mu_weighted(W | F, Uedges, al, A, memo))
        Ws.add(W)
    I = A.ideal([e for e in eqs if e != 0])
    for W in Ws:
        f = dvg.mu_weighted(W | F, Uedges, al, A, memo)
        if f.is_constant():
            if f == 0: return False
            continue
        I = I.saturation(A.ideal([f]))[0]
        if I.is_one(): return False
    return not I.is_one()

if __name__ == '__main__':
    n = 3
    tests = {'C3': [(0,1),(1,2),(2,0)], 'tri 0<->1,0<->2,1->2': [(0,1),(1,0),(0,2),(2,0),(1,2)],
             'tri 0<->1,1->2,2->0': [(0,1),(1,0),(1,2),(2,0)]}
    for name, arcs in tests.items():
        for nf in (1, 2):
            N = n + nf
            pairs = [(a, b) for a in range(N) for b in range(a + 1, N)]
            found = []
            for mask in range(1 << len(pairs)):
                U = frozenset(pairs[i] for i in range(len(pairs)) if mask >> i & 1)
                if has_env_potential(n, arcs, nf, U):
                    found.append(sorted(U))
            print(name, '|F| =', nf, 'U with environment potential:', len(found), found[:6], flush=True)
