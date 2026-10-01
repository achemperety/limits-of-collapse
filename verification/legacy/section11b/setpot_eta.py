"""Set potentials with self-energies, decided exactly over Q(z).

setpot_eta(n, arcs, eta) plays the game with 1/R(W, v) = z - eta[v] - sum R(W - v, t) and propagates
M(W - v) = R(W, v) M(W) from M(V) = 1 along the reachable positions; it returns True exactly when every
reachable unvisited set receives a single value, that is, when Gamma(P) depends only on the vertex set of P.
"""
from functools import lru_cache
import dvg

Kz, z = dvg.Kz, dvg.Kz(dvg.z)


def setpot_eta(n, arcs, eta):
    out = dvg.out_nbrs(n, arcs)
    eta = {v: Kz(eta.get(v, 0)) for v in range(n)}

    @lru_cache(maxsize=None)
    def R(W, v):
        s = eta[v]
        for t in out[v]:
            if t in W and t != v:
                s += R(W - {v}, t)
        return 1 / (z - s)
    P = dvg.reachable_positions(n, arcs)
    M = {frozenset(range(n)): Kz(1)}
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


def cycles_digraph(cycles):
    """Union of directed cycles given as vertex lists (labels 0..n-1)."""
    arcs = set()
    for c in cycles:
        for i in range(len(c)):
            arcs.add((c[i], c[(i + 1) % len(c)]))
    n = 1 + max(max(c) for c in cycles)
    return n, sorted(arcs)
