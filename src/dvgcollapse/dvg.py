"""Directed vertex geography (DVG): digraphs, positions, outcomes and exact resolvents.

Conventions (as in the paper):
  * A digraph D has vertex set {0, ..., n-1} and a set of arcs (u, v), u != v.
  * A position is a pair (W, v): W is the frozenset of unvisited vertices, v in W is the token.
    The options of (W, v) are (W - {v}, w) for out-neighbours w of v in W.
  * The player to move at a position with no option loses (normal play).
  * The root resolvent R(W, v) is the resolvent at the root of the game tree; it satisfies
        1/R(W, v) = z - eta(v) - sum_{w} R(W - {v}, w)
    (Theorem 5.2 and the self-energy convention of Section 9), where eta is an optional
    self-energy attached to each vertex (default 0).

Everything here is exact: outcomes are Booleans, resolvents are elements of the rational
function field Q(z) (sympy.polys.fields), or rationals when z is given as a Fraction.
"""
from __future__ import annotations

import itertools
from fractions import Fraction
from functools import lru_cache
from typing import Dict, FrozenSet, Iterable, List, Optional, Sequence, Set, Tuple

from sympy import QQ, field

# The rational function field Q(z) used throughout.
QZ, Z = field('z', QQ)

Position = Tuple[FrozenSet[int], int]


class Digraph:
    """A finite simple digraph on {0, ..., n-1}."""

    def __init__(self, n: int, arcs: Iterable[Tuple[int, int]]):
        self.n = int(n)
        self.arcs: Set[Tuple[int, int]] = set()
        for a, b in arcs:
            if a == b:
                raise ValueError('loops are not allowed')
            if not (0 <= a < n and 0 <= b < n):
                raise ValueError('arc (%r, %r) out of range' % (a, b))
            self.arcs.add((int(a), int(b)))
        self.out: List[List[int]] = [sorted(b for (a, b) in self.arcs if a == v) for v in range(n)]
        self.inn: List[List[int]] = [sorted(a for (a, b) in self.arcs if b == v) for v in range(n)]

    # ------------------------------------------------------------------ structure
    def has_arc(self, a: int, b: int) -> bool:
        return (a, b) in self.arcs

    def reach(self) -> List[Set[int]]:
        """reach()[v] = set of vertices reachable from v (including v)."""
        res = []
        for v in range(self.n):
            seen = {v}
            stack = [v]
            while stack:
                x = stack.pop()
                for y in self.out[x]:
                    if y not in seen:
                        seen.add(y)
                        stack.append(y)
            res.append(seen)
        return res

    def strong_components(self) -> List[FrozenSet[int]]:
        """Strong components in reverse topological order of the condensation
        (sinks first), computed with Tarjan's algorithm."""
        index: Dict[int, int] = {}
        low: Dict[int, int] = {}
        on: Set[int] = set()
        stack: List[int] = []
        comps: List[FrozenSet[int]] = []
        counter = [0]

        def visit(v: int) -> None:
            # iterative Tarjan to avoid recursion limits
            work = [(v, 0)]
            index[v] = low[v] = counter[0]
            counter[0] += 1
            stack.append(v)
            on.add(v)
            while work:
                x, i = work[-1]
                if i < len(self.out[x]):
                    work[-1] = (x, i + 1)
                    y = self.out[x][i]
                    if y not in index:
                        index[y] = low[y] = counter[0]
                        counter[0] += 1
                        stack.append(y)
                        on.add(y)
                        work.append((y, 0))
                    elif y in on:
                        low[x] = min(low[x], index[y])
                else:
                    work.pop()
                    if work:
                        p = work[-1][0]
                        low[p] = min(low[p], low[x])
                    if low[x] == index[x]:
                        comp = set()
                        while True:
                            w = stack.pop()
                            on.discard(w)
                            comp.add(w)
                            if w == x:
                                break
                        comps.append(frozenset(comp))

        for v in range(self.n):
            if v not in index:
                visit(v)
        return comps  # Tarjan emits components in reverse topological order

    def component_of(self) -> Dict[int, FrozenSet[int]]:
        return {v: C for C in self.strong_components() for v in C}

    def is_scc_symmetric(self) -> bool:
        """Every arc lying on a directed cycle has its reverse arc (Section 6.5)."""
        comp = self.component_of()
        return all((b, a) in self.arcs for (a, b) in self.arcs if comp[a] is comp[b] or comp[a] == comp[b])

    def one_way_arcs_in_components(self) -> List[Tuple[int, int]]:
        comp = self.component_of()
        return sorted((a, b) for (a, b) in self.arcs if comp[a] == comp[b] and (b, a) not in self.arcs)

    def symmetric_part(self) -> Set[Tuple[int, int]]:
        """Edges {a, b} (as pairs a < b) joined in both directions."""
        return {(min(a, b), max(a, b)) for (a, b) in self.arcs if (b, a) in self.arcs}

    def induced(self, S: Iterable[int]) -> Tuple['Digraph', List[int]]:
        """Induced subdigraph on S, relabelled 0..|S|-1; returns (digraph, old labels)."""
        S = sorted(S)
        pos = {v: i for i, v in enumerate(S)}
        return Digraph(len(S), [(pos[a], pos[b]) for (a, b) in self.arcs if a in pos and b in pos]), S

    # ------------------------------------------------------------------ game
    def options(self, W: FrozenSet[int], v: int) -> List[int]:
        return [w for w in self.out[v] if w in W and w != v]

    def reachable_positions(self, starts: Optional[Iterable[int]] = None) -> Set[Position]:
        """All positions reachable from the starts (V, s) (all s by default)."""
        V = frozenset(range(self.n))
        seen: Set[Position] = set()
        todo = [(V, s) for s in (range(self.n) if starts is None else starts)]
        while todo:
            W, v = todo.pop()
            if (W, v) in seen:
                continue
            seen.add((W, v))
            for w in self.options(W, v):
                todo.append((W - {v}, w))
        return seen

    def outcome_function(self):
        """Return win(W, v): True iff the player to move at (W, v) wins (N-position)."""
        out = self.out

        @lru_cache(maxsize=None)
        def win(W: FrozenSet[int], v: int) -> bool:
            W2 = W - {v}
            return any(not win(W2, w) for w in out[v] if w in W2)

        return win

    def wins_from(self, v: int, W: Optional[Iterable[int]] = None) -> bool:
        W = frozenset(range(self.n)) if W is None else frozenset(W)
        return self.outcome_function()(W, v)

    def resolvent_function(self, eta: Optional[Dict[int, object]] = None, zval=None):
        """Return R(W, v), the exact root resolvent.

        With zval=None the values lie in Q(z) (sympy FracElement); with zval a Fraction
        the values are Fractions (z specialised).  eta maps vertices to self-energies in
        the same field (default 0)."""
        out = self.out
        if zval is None:
            zz = Z
            zero = QZ(0)
        else:
            zz = Fraction(zval)
            zero = Fraction(0)
        et = {v: (eta.get(v, zero) if eta else zero) for v in range(self.n)}

        @lru_cache(maxsize=None)
        def R(W: FrozenSet[int], v: int):
            W2 = W - {v}
            s = zero
            for w in out[v]:
                if w in W2:
                    s = s + R(W2, w)
            return 1 / (zz - et[v] - s)

        return R

    def gamma(self, path: Sequence[int], R=None):
        """Gamma(P): product of R(V - {t_0..t_{i-1}}, t_i) along a directed path P started at V."""
        if R is None:
            R = self.resolvent_function()
        W = frozenset(range(self.n))
        g = None
        for i, t in enumerate(path):
            if i > 0 and not self.has_arc(path[i - 1], t):
                raise ValueError('not a directed path')
            if t not in W:
                raise ValueError('path revisits a vertex')
            r = R(W, t)
            g = r if g is None else g * r
            W = W - {t}
        return g

    def back_arcs(self, path: Sequence[int]) -> int:
        """back(P): number of arcs t_j -> t_i with i < j (Lemma 9.7; Lemma 10.A of the addendum)."""
        return sum(1 for i in range(len(path)) for j in range(i + 1, len(path)) if self.has_arc(path[j], path[i]))


def order_at_zero(f) -> int:
    """Order of vanishing at z = 0 of a nonzero element of Q(z) (negative for a pole)."""
    num, den = f.numer, f.denom

    def ordp(p) -> int:
        # p is a PolyElement in z; lowest exponent with nonzero coefficient
        terms = p.terms()
        return min(m[0] for m, c in terms)

    return ordp(num) - ordp(den)


def all_digraphs(n: int):
    """All digraphs on n labelled vertices up to isomorphism (brute force; n <= 4)."""
    pairs = [(a, b) for a in range(n) for b in range(n) if a != b]
    seen = set()
    for mask in range(1 << len(pairs)):
        arcs = [pairs[i] for i in range(len(pairs)) if mask >> i & 1]
        best = None
        for p in itertools.permutations(range(n)):
            key = tuple(sorted((p[a], p[b]) for a, b in arcs))
            if best is None or key < best:
                best = key
        if best in seen:
            continue
        seen.add(best)
        yield Digraph(n, best)


def random_digraph(n: int, p: float, rng) -> Digraph:
    return Digraph(n, [(a, b) for a in range(n) for b in range(n) if a != b and rng.random() < p])
