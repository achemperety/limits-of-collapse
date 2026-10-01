"""The common-tail ("star-tail") evaluator: Theorem 12.2 and Corollaries 12.3-12.4
(Theorems 10.11, 11.1 and Corollary 11.1' of the addendum).

Game Gamma(H; A).  H is an undirected graph; A is a set of one-way arcs u -> w between
non-adjacent vertices of H.  The token moves along an edge of H to an unvisited vertex, or
along an arc of A whose head is unvisited; the vertex it leaves is deleted; a player who
cannot move loses.  This is the game played inside a strong component of a digraph once
exits have been replaced by pendant leaves (Section 12.1, "Reduction to one component").

Star rule (Theorem 12.2).  Let every arc of A leave the same tail u.  At a position with
current graph K (the unvisited vertices, with the edges of H among them), u unvisited, token
on s, let W' be the set of live heads other than s and K + uW' the graph K with the edges
uw (w in W') added.  K is *split* if nu(K + uW') = nu(K), *joined* otherwise.  The player to
move wins iff
    1. s = u and u is essential in K + uW'; or
    2. s != u, s is essential in K, and (K - s is split or u is inessential in K); or
    3. s != u, s is inessential in K, K is joined and nu(K - s - u) < nu(K).
With no live arc the rule is Theorem 6.2: the player to move wins iff s is essential in K.

Every quantity is a matching number of K minus at most three key points, because
    nu(K' + uW') = max(nu(K'), max_{w in W'} nu(K' - u - w) + 1)             (equation (12.1)),
so each position costs O(|W'|) maximum-matching computations.

The module also provides
    * gamma_outcome_bruteforce   exhaustive evaluation of Gamma(H; A) (for testing),
    * solve_dvg                  DVG on a digraph whose strong components are symmetric or have
                                 one-way arcs with a common tail (Corollary 12.4), for any position.
"""
from __future__ import annotations

from functools import lru_cache
from typing import Dict, FrozenSet, Iterable, List, Optional, Sequence, Set, Tuple

from .matching import nu as _nu
from .dvg import Digraph

Edge = Tuple[int, int]
Arc = Tuple[int, int]


class NotStarShaped(ValueError):
    """Raised when the one-way arcs of a strong component do not share a tail."""


# ---------------------------------------------------------------------------- matching helpers
class _Matcher:
    """Memoised matching numbers of induced subgraphs of one fixed graph."""

    def __init__(self, vertices: Iterable[int], edges: Iterable[Edge]):
        self.V = frozenset(vertices)
        self.E = [(a, b) for (a, b) in edges if a in self.V and b in self.V and a != b]
        self._memo: Dict[FrozenSet[int], int] = {}

    def nu(self, S: Iterable[int]) -> int:
        S = frozenset(S)
        if S not in self._memo:
            self._memo[S] = _nu(S, self.E)
        return self._memo[S]

    def nu_plus_star(self, S: FrozenSet[int], u: int, heads: Iterable[int]) -> int:
        """nu(K + uW') for K = G[S], u in S, W' = heads (subset of S - {u})."""
        best = self.nu(S)
        for w in heads:
            best = max(best, self.nu(S - {u, w}) + 1)
        return best


def star_rule(K: Iterable[int], edges: Iterable[Edge], u: Optional[int], heads: Iterable[int],
              s: int, matcher: Optional[_Matcher] = None) -> bool:
    """Outcome (True = player to move wins) of the position of Gamma(H; A) with current vertex
    set K (edges of H among K are taken from `edges`), token on s, common tail u (None if the
    tail is visited or there are no arcs) and live heads `heads` (heads in K; s may be among them).

    Implements Theorem 12.2 (Corollary 12.3 when there is one head)."""
    K = frozenset(K)
    if s not in K:
        raise ValueError('token must lie in K')
    m = matcher or _Matcher(K, edges)
    Wp = frozenset(h for h in heads if h in K and h != s)
    if u is None or u not in K or not Wp:
        # plain undirected vertex geography (Theorem 6.2)
        return m.nu(K - {s}) == m.nu(K) - 1
    if s == u:
        # item 1: u essential in K + uW', i.e. nu((K + uW') - u) = nu(K + uW') - 1,
        # and (K + uW') - u = K - u because every added edge contains u.
        return m.nu(K - {u}) == m.nu_plus_star(K, u, Wp) - 1
    nK = m.nu(K)
    s_essential = m.nu(K - {s}) == nK - 1
    if s_essential:
        Ks = K - {s}
        split_Ks = m.nu_plus_star(Ks, u, Wp) == m.nu(Ks)
        u_inessential = m.nu(K - {u}) == nK
        return split_Ks or u_inessential
    joined = m.nu_plus_star(K, u, Wp) > nK
    return joined and m.nu(K - {s, u}) < nK


def _check_star(arcs: Iterable[Arc]) -> Optional[int]:
    tails = {a for (a, b) in arcs}
    if len(tails) > 1:
        raise NotStarShaped('one-way arcs with tails %s' % sorted(tails))
    return next(iter(tails)) if tails else None


def gamma_outcome(H_vertices: Iterable[int], H_edges: Iterable[Edge], arcs: Iterable[Arc],
                  s: int, K: Optional[Iterable[int]] = None) -> bool:
    """Outcome of Gamma(H; A) at the position (K, s) (K defaults to all vertices of H) by the
    star rule; the arcs must share a tail."""
    arcs = list(arcs)
    u = _check_star(arcs)
    Kset = frozenset(H_vertices if K is None else K)
    heads = [w for (a, w) in arcs]
    return star_rule(Kset, H_edges, u, heads, s)


def gamma_outcome_bruteforce(H_vertices: Iterable[int], H_edges: Iterable[Edge], arcs: Iterable[Arc],
                             s: int, K: Optional[Iterable[int]] = None) -> bool:
    """Exhaustive game-tree evaluation of Gamma(H; A) at (K, s) (any arc geometry)."""
    adj: Dict[int, Set[int]] = {}
    for a, b in H_edges:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    arc_out: Dict[int, Set[int]] = {}
    for a, b in arcs:
        arc_out.setdefault(a, set()).add(b)
    K0 = frozenset(H_vertices if K is None else K)

    @lru_cache(maxsize=None)
    def win(W: FrozenSet[int], v: int) -> bool:
        W2 = W - {v}
        moves = (adj.get(v, set()) | arc_out.get(v, set())) & W2
        return any(not win(W2, w) for w in moves)

    return win(K0, s)


# ---------------------------------------------------------------------------- DVG reduction
def _component_game(D: Digraph, C: FrozenSet[int], fresh_loses: Dict[int, bool]):
    """Build Gamma(H; A) for the strong component C: H = symmetric pairs of C plus one pendant
    leaf at every vertex with an exit to a P-position (fresh entry that loses for the mover),
    A = one-way arcs inside C.  Returns (H_vertices, H_edges, arcs)."""
    Hv = set(C)
    He: List[Edge] = []
    for (a, b) in D.arcs:
        if a in C and b in C and a < b and (b, a) in D.arcs:
            He.append((a, b))
    arcs = [(a, b) for (a, b) in D.arcs if a in C and b in C and (b, a) not in D.arcs]
    leaf = D.n  # fresh labels for pendant leaves
    for x in sorted(C):
        if any(y not in C and fresh_loses[y] for y in D.out[x]):
            Hv.add(leaf)
            He.append((x, leaf))
            leaf += 1
    return Hv, He, arcs


def fresh_entry_outcomes(D: Digraph) -> Dict[int, bool]:
    """win[y] = True iff the player to move wins DVG on D started at y (all vertices unvisited),
    computed component by component in reverse topological order (Corollaries 6.8, 10.11', 11.1').
    Raises NotStarShaped if some component's one-way arcs have more than one tail."""
    win: Dict[int, bool] = {}
    loses: Dict[int, bool] = {}
    for C in D.strong_components():  # sinks first
        Hv, He, arcs = _component_game(D, C, loses)
        u = _check_star(arcs)
        heads = [w for (a, w) in arcs]
        m = _Matcher(Hv, He)
        for y in C:
            win[y] = star_rule(frozenset(Hv), He, u, heads, y, matcher=m)
            loses[y] = not win[y]
    return win


def solve_dvg(D: Digraph, v: int, W: Optional[Iterable[int]] = None) -> bool:
    """Outcome of the DVG position (W, v) of D (W defaults to all vertices): DVG from (W, v)
    is DVG on D[W] started at v.  Requires the star condition in every strong component of D[W]."""
    W = frozenset(range(D.n)) if W is None else frozenset(W)
    if v not in W:
        raise ValueError('token must be unvisited')
    sub, labels = D.induced(W)
    pos = {x: i for i, x in enumerate(labels)}
    return fresh_entry_outcomes(sub)[pos[v]]


def key_local_data(K: Iterable[int], edges: Iterable[Edge], keys: Sequence[int]) -> Tuple:
    """Key-local data (Definition 12.1): coincidence pattern of the key points, the deficiencies
    nu(K) - nu(K - X) for all sets X of distinct key points, and the adjacency among key points."""
    K = frozenset(K)
    m = _Matcher(K, edges)
    E = {frozenset(e) for e in edges}
    distinct = sorted(set(keys), key=list(keys).index)
    pattern = tuple(distinct.index(k) for k in keys)
    n0 = m.nu(K)
    defs = []
    for r in range(1, len(distinct) + 1):
        from itertools import combinations
        for X in combinations(range(len(distinct)), r):
            defs.append(n0 - m.nu(K - {distinct[i] for i in X}))
    adj = tuple(int(frozenset((distinct[i], distinct[j])) in E)
                for i in range(len(distinct)) for j in range(i + 1, len(distinct)))
    return pattern, tuple(defs), adj
