"""Ladders, intrinsic degree, canonical subgames, the memory comb and component-state (circuit)
potentials (Section 14 of the revision).

  * ladder(m)                        the SCC-symmetric ladder L_m of Theorem 14.1 (pathwidth 2)
  * closed_ladder(m)                 L_m plus the arc b_{m-1} -> t_0 (strongly connected, not symmetric)
  * canonical_subgame(D, W, v)       (reach of v in D[W], v): the state on which the recursion closes (Lemma 14.4)
  * comb(m)                          the treewidth-3 comb of Proposition 14.5 with prime tails
  * resolvent_degree(D, s)           degree of the denominator of R(V, s) in lowest terms
  * ComponentStatePotential(D)       exact resolvents of every reachable position from the
                                     component states (C, S, v), S subset of C; the recursion is an
                                     arithmetic circuit with sum_C |C| 2^{|C|} gates
"""
from __future__ import annotations

from functools import lru_cache
from typing import Dict, FrozenSet, List, Tuple

from .dvg import Digraph, QZ, Z


def ladder(m: int) -> Digraph:
    """t_i = i, b_i = m + i; arcs t_i -> t_{i+1}, b_i -> b_{i+1}; rungs t_i <-> b_i."""
    arcs = []
    for i in range(m - 1):
        arcs += [(i, i + 1), (m + i, m + i + 1)]
    for i in range(m):
        arcs += [(i, m + i), (m + i, i)]
    return Digraph(2 * m, arcs)


def resolvent_degree(D: Digraph, s: int) -> int:
    R = D.resolvent_function()
    r = R(frozenset(range(D.n)), s)
    return r.denom.degree()


class ComponentStatePotential:
    """Exact resolvents through component states.

    Fresh-entry resolvents rho(y) = R(V, y) are computed component by component in reverse
    topological order; inside a component C with exit self-energies eta_x = sum_{exits y} rho(y),
        R_C(S, v) = 1 / (z - eta_v - sum_{w in N^+(v) cap (S - v)} R_C(S - v, w)),   S subset of C.
    By the entry lemma (Lemma 6.5), every reachable position (W, v) with v in C satisfies
        R(W, v) = R_C(W cap C, v).
    The number of states is sum_C |C| 2^{|C|}, polynomial when every strong component has
    O(log n) vertices, although the number of reachable positions can be exponential (ladders)."""

    def __init__(self, D: Digraph):
        self.D = D
        self.comps = D.strong_components()          # sinks first
        self.comp = {v: C for C in self.comps for v in C}
        self.rho: Dict[int, object] = {}
        self._tables: Dict[FrozenSet[int], object] = {}
        for C in self.comps:
            eta = {x: sum((self.rho[y] for y in D.out[x] if y not in C), QZ(0)) for x in C}
            RC = self._component_resolvent(C, eta)
            self._tables[C] = RC
            full = frozenset(C)
            for y in C:
                self.rho[y] = RC(full, y)

    def _component_resolvent(self, C: FrozenSet[int], eta):
        out = self.D.out

        @lru_cache(maxsize=None)
        def RC(S: FrozenSet[int], v: int):
            S2 = S - {v}
            acc = QZ(0)
            for w in out[v]:
                if w in S2:
                    acc = acc + RC(S2, w)
            return 1 / (Z - eta[v] - acc)
        return RC

    def resolvent(self, W, v):
        C = self.comp[v]
        return self._tables[C](frozenset(W) & C, v)

    def number_of_states(self) -> int:
        return sum(len(C) * (1 << len(C)) for C in self.comps)


def closed_ladder(m: int) -> Digraph:
    """L_m plus the arc b_{m-1} -> t_0: strongly connected, one-way rails, pathwidth at most 3."""
    L = ladder(m)
    return Digraph(2 * m, set(L.arcs) | {(2 * m - 1, 0)})


def canonical_subgame(D: Digraph, W, v: int) -> Tuple[FrozenSet[int], int]:
    """The set of vertices reachable from v inside D[W], with the token v (Lemma 14.4)."""
    W = frozenset(W)
    seen = {v}
    stack = [v]
    while stack:
        x = stack.pop()
        for y in D.out[x]:
            if y in W and y not in seen:
                seen.add(y)
                stack.append(y)
    return frozenset(seen), v


def _odd_primes(k: int) -> List[int]:
    out, c = [], 3
    while len(out) < k:
        if all(c % p for p in range(3, int(c ** 0.5) + 1, 2)):
            out.append(c)
        c += 2
    return out


def comb(m: int, primes: List[int] = None):
    """The comb Q_m of Proposition 14.5.

    Spine p_0..p_m; for each i < m two vertices a_i, b_i with arcs p_i -> a_i -> p_{i+1} and
    p_i -> b_i -> p_{i+1}; a collector q with p_m -> q and q -> a_i, q -> b_i; and at every
    x in {a_i, b_i} a private directed tail, so that the directed path from x through its tail has
    L(x) vertices with L(x) + 1 prime, all 2m primes distinct.  Returns (D, index) where index maps
    names ('p', i), ('a', i), ('b', i), 'q' to vertex numbers."""
    primes = primes or _odd_primes(2 * m)
    idx: Dict[object, int] = {}

    def v(name):
        if name not in idx:
            idx[name] = len(idx)
        return idx[name]

    arcs = []
    for i in range(m):
        arcs += [(v(('p', i)), v(('a', i))), (v(('a', i)), v(('p', i + 1))),
                 (v(('p', i)), v(('b', i))), (v(('b', i)), v(('p', i + 1)))]
        arcs += [(v('q'), v(('a', i))), (v('q'), v(('b', i)))]
    arcs.append((v(('p', m)), v('q')))
    for i in range(m):
        for j, lab in enumerate(('a', 'b')):
            L = primes[2 * i + j] - 1           # vertices on the path from x through its tail
            prev = v((lab, i))
            for k in range(L - 1):
                t = v((lab + 't', i, k))
                arcs.append((prev, t))
                prev = t
    return Digraph(len(idx), arcs), idx


def resolvents_at_collector(m: int):
    """The set of resolvents R(W, q) over the positions at q reachable from (V, p_0)."""
    D, idx = comb(m)
    R = D.resolvent_function()
    P = D.reachable_positions(starts=[idx[('p', 0)]])
    return {R(W, idx['q']) for (W, v) in P if v == idx['q']}
