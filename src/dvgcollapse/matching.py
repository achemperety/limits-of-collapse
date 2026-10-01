"""Matchings: matching numbers, Gallai-Edmonds classes and weighted matching polynomials.

  * nu(vertices, edges)             maximum matching number (Edmonds' blossom algorithm via networkx)
  * essential(v, vertices, edges)   v is covered by every maximum matching  (Theorem 6.2)
  * gallai_edmonds(vertices, edges) the classes D, A, C of the Gallai-Edmonds decomposition
  * mu_weighted(S, edges, a)        mu_a(U[S]) = sum over matchings M of U[S] of (-1)^|M| prod_{x not in V(M)} a_x
"""
from __future__ import annotations

from typing import Dict, FrozenSet, Iterable, Set, Tuple

import networkx as nx

Edge = Tuple[int, int]


def _graph(vertices: Iterable[int], edges: Iterable[Edge]) -> nx.Graph:
    V = set(vertices)
    G = nx.Graph()
    G.add_nodes_from(V)
    G.add_edges_from((a, b) for (a, b) in edges if a in V and b in V and a != b)
    return G


def nu(vertices: Iterable[int], edges: Iterable[Edge]) -> int:
    """Maximum matching number of the graph (vertices, edges restricted to vertices)."""
    G = _graph(vertices, edges)
    if G.number_of_edges() == 0:
        return 0
    return len(nx.max_weight_matching(G, maxcardinality=True))


def essential(v: int, vertices: Iterable[int], edges: Iterable[Edge]) -> bool:
    """True iff v is covered by every maximum matching, i.e. nu(G - v) = nu(G) - 1."""
    V = set(vertices)
    E = list(edges)
    return nu(V - {v}, E) == nu(V, E) - 1


def gallai_edmonds(vertices: Iterable[int], edges: Iterable[Edge]) -> Dict[str, Set[int]]:
    """Gallai-Edmonds decomposition: D = vertices missed by some maximum matching,
    A = neighbours of D outside D, C = the rest."""
    V = set(vertices)
    E = list(edges)
    n0 = nu(V, E)
    D = {v for v in V if nu(V - {v}, E) == n0}
    adj = {v: set() for v in V}
    for a, b in E:
        if a in V and b in V:
            adj[a].add(b)
            adj[b].add(a)
    A = {v for v in V - D if adj[v] & D}
    return {'D': D, 'A': A, 'C': V - D - A}


def mu_weighted(S: Iterable[int], edges: Iterable[Edge], a: Dict[int, object], one=1, memo=None):
    """Weighted matching polynomial mu_a(U[S]) by the vertex recurrence
        mu(H) = a_v mu(H - v) - sum_{w ~ v} mu(H - v - w).
    The values a_v may lie in any commutative ring (Fractions, sympy field elements, ...);
    `one` is the unit of that ring.  Edges are pairs (i, j); orientation is ignored."""
    S = frozenset(S)
    if memo is None:
        memo = {}
    adj: Dict[int, Set[int]] = {}
    for (i, j) in edges:
        adj.setdefault(i, set()).add(j)
        adj.setdefault(j, set()).add(i)

    def rec(T: FrozenSet[int]):
        if T in memo:
            return memo[T]
        if not T:
            return one
        v = min(T)
        rest = T - {v}
        val = a[v] * rec(rest)
        for w in adj.get(v, ()):
            if w in rest:
                val = val - rec(rest - {w})
        memo[T] = val
        return val

    return rec(S)


def matching_polynomial(vertices: Iterable[int], edges: Iterable[Edge]):
    """Ordinary matching polynomial mu(G, z) in Q(z) (all weights equal to z)."""
    from .dvg import QZ, Z
    V = list(vertices)
    return mu_weighted(V, edges, {v: Z for v in V}, one=QZ(1))
