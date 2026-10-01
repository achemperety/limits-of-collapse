"""Potentials for directed vertex geography.

  * scc_symmetric_potential(D)   the local weighted matching potential of Theorem 6.6
  * verify_potential(D, U, a)    exact check of mu_a(U[W - v]) = R(W, v) mu_a(U[W]), mu_a(U[W]) != 0
  * set_potential(D, eta)        M with M(W - v) = R(W, v) M(W) on reachable positions, or None
  * valuation_potential(D)       Lambda with Lambda(W - v) - Lambda(W) = ord_0 R(W, v), or None
"""
from __future__ import annotations

from typing import Dict, FrozenSet, Iterable, Optional, Set, Tuple

from .dvg import Digraph, QZ, Z, order_at_zero
from .matching import mu_weighted

Edge = Tuple[int, int]


def scc_symmetric_potential(D: Digraph):
    """Theorem 6.6: for SCC-symmetric D, U = symmetric part and a_x = z - eta_x with
    eta_x = sum of the fresh-entry resolvents R(V, y) over exits y of x (out-neighbours in
    other strong components).  Returns (U, a) with a[x] in Q(z)."""
    if not D.is_scc_symmetric():
        raise ValueError('Theorem 6.6 needs an SCC-symmetric digraph')
    comp = D.component_of()
    R = D.resolvent_function()
    V = frozenset(range(D.n))
    U = D.symmetric_part()
    a = {}
    for x in range(D.n):
        eta = QZ(0)
        for y in D.out[x]:
            if comp[y] != comp[x]:
                eta = eta + R(V, y)
        a[x] = Z - eta
    return U, a


def verify_potential(D: Digraph, U: Iterable[Edge], a: Dict[int, object], positions=None) -> bool:
    """Exact verification of a local weighted matching potential (U, a) over Q(z)."""
    U = list(U)
    R = D.resolvent_function()
    memo: Dict[FrozenSet[int], object] = {}
    for (W, v) in (positions or D.reachable_positions()):
        den = mu_weighted(W, U, a, one=QZ(1), memo=memo)
        if den == 0:
            return False
        if mu_weighted(W - {v}, U, a, one=QZ(1), memo=memo) != R(W, v) * den:
            return False
    return True


def set_potential(D: Digraph, eta: Optional[Dict[int, object]] = None):
    """Propagate M(V) = 1 along M(W - v) = R(W, v) M(W); return the dict M or None if some
    reachable unvisited set receives two different values (then no set potential exists over
    any field containing Q(z))."""
    R = D.resolvent_function(eta=eta)
    V = frozenset(range(D.n))
    M: Dict[FrozenSet[int], object] = {V: QZ(1)}
    for (W, v) in sorted(D.reachable_positions(), key=lambda p: -len(p[0])):
        if W not in M:
            continue
        val = M[W] * R(W, v)
        W2 = W - {v}
        if W2 in M:
            if M[W2] != val:
                return None
        else:
            M[W2] = val
    return M


def valuation_potential(D: Digraph):
    """Integer Lambda with Lambda(W - v) - Lambda(W) = +1 at N-positions and -1 at P-positions
    (Theorem 5.2), or None."""
    win = D.outcome_function()
    V = frozenset(range(D.n))
    L: Dict[FrozenSet[int], int] = {V: 0}
    for (W, v) in sorted(D.reachable_positions(), key=lambda p: -len(p[0])):
        if W not in L:
            continue
        val = L[W] + (1 if win(W, v) else -1)
        W2 = W - {v}
        if W2 in L:
            if L[W2] != val:
                return None
        else:
            L[W2] = val
    return L


def check_outcomes_by_valuation(D: Digraph) -> bool:
    """Theorem 5.2 on every reachable position: ord_0 R(W, v) = +1 iff (W, v) is an N-position."""
    R = D.resolvent_function()
    win = D.outcome_function()
    return all((order_at_zero(R(W, v)) == 1) == win(W, v) and order_at_zero(R(W, v)) in (1, -1)
               for (W, v) in D.reachable_positions())
