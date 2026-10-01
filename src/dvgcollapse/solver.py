"""A unified solver for directed vertex geography (DVG), dispatching to the exact methods of the paper.

    outcome(D, v, W=None)  ->  (wins, method)

Methods, tried in this order (every one is exact):
  'common-tail'       every strong component of D[W] has one-way arcs with at most one tail
                      (Corollary 6.8 for symmetric components, Theorem 12.2 / Corollary 12.4 otherwise):
                      polynomial time, at most six maximum-matching computations per position;
  'component-states'  every non-symmetric strong component has at most `max_component` vertices
                      (Corollary 14.7): time 2^k poly(n);
  'exhaustive'        memoized search over unvisited sets (exponential; always correct).

Command line:
    python -m dvgcollapse.solver --n 4 --arcs "0-1 1-2 2-0 2-3" --start 0 [--unvisited "0 1 2 3"]
"""
from __future__ import annotations

import argparse
from functools import lru_cache
from typing import Dict, FrozenSet, Iterable, Optional, Tuple

from .dvg import Digraph
from .matching import essential
from .star_tail import NotStarShaped, solve_dvg


def _component_state_outcome(D: Digraph, v: int, max_component: int) -> bool:
    """Outcome of the fresh start (V, v) by the component-state recursion (Corollary 14.7)."""
    comps = D.strong_components()             # sinks first
    comp = {x: C for C in comps for x in C}
    fresh: Dict[int, bool] = {}               # outcome of a fresh entry at each vertex
    for C in comps:
        one_way = any(not D.has_arc(b, a) for a in C for b in D.out[a] if b in C)
        if one_way and len(C) > max_component:
            raise NotStarShaped('component too large for the component-state recursion')
        exits_win = {x: any((y not in C) and not fresh[y] for y in D.out[x]) for x in C}
        if not one_way:
            # symmetric component: undirected geography on G(C) plus a pendant leaf at every
            # vertex with an exit to a P-position (Corollary 6.8); the mover wins iff essential
            edges = [(a, b) for a in C for b in D.out[a] if b in C and a < b]
            verts = list(C)
            leaf = max(D.n, max(C) + 1)
            for x in C:
                if exits_win[x]:
                    edges.append((x, leaf))
                    verts.append(leaf)
                    leaf += 1
            for x in C:
                fresh[x] = essential(x, verts, edges)
            continue

        @lru_cache(maxsize=None)
        def win(S: FrozenSet[int], x: int, _C=C, _ew=exits_win) -> bool:
            if _ew[x]:
                return True
            S2 = S - {x}
            return any(not win(S2, y) for y in D.out[x] if y in S2)

        full = frozenset(C)
        for x in C:
            fresh[x] = win(full, x)
    return fresh[v]


def outcome(D: Digraph, v: int, W: Optional[Iterable[int]] = None,
            max_component: int = 16) -> Tuple[bool, str]:
    """Return (True if the player to move at (W, v) wins, name of the method used)."""
    W = frozenset(range(D.n)) if W is None else frozenset(W)
    if v not in W:
        raise ValueError('the token must be on an unvisited vertex')
    try:
        return solve_dvg(D, v, W), 'common-tail'
    except NotStarShaped:
        pass
    sub, labels = D.induced(W)
    pos = {x: i for i, x in enumerate(labels)}
    try:
        return _component_state_outcome(sub, pos[v], max_component), 'component-states'
    except NotStarShaped:
        pass
    return D.outcome_function()(W, v), 'exhaustive'


def _parse_arcs(text: str):
    arcs = []
    for tok in text.replace(',', ' ').split():
        a, b = tok.split('-')
        arcs.append((int(a), int(b)))
    return arcs


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description='Exact outcome of a directed vertex geography position.')
    ap.add_argument('--n', type=int, required=True, help='number of vertices (0..n-1)')
    ap.add_argument('--arcs', required=True, help='arcs as "a-b a-b ..."')
    ap.add_argument('--start', type=int, required=True, help='vertex holding the token')
    ap.add_argument('--unvisited', default=None, help='unvisited vertices (default: all)')
    ap.add_argument('--max-component', type=int, default=16)
    args = ap.parse_args(argv)
    D = Digraph(args.n, _parse_arcs(args.arcs))
    W = None if args.unvisited is None else [int(x) for x in args.unvisited.replace(',', ' ').split()]
    wins, method = outcome(D, args.start, W, args.max_component)
    print('%s (method: %s)' % ('first player wins' if wins else 'first player loses', method))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
