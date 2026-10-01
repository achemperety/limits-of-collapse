"""Tests for the common-tail evaluator (Theorem 12.2, Corollaries 12.3-12.4)."""
import itertools
import random

import pytest

from dvgcollapse.dvg import Digraph
from dvgcollapse.star_tail import (NotStarShaped, gamma_outcome, gamma_outcome_bruteforce, key_local_data,
                                   solve_dvg, star_rule)


def _edges(spec):
    return [tuple(map(int, e.split('-'))) for e in spec.split()]


def _random_instance(rng, n, heads_max):
    V = list(range(n))
    E = [(a, b) for a in V for b in V if a < b and rng.random() < 0.35]
    Es = {frozenset(e) for e in E}
    u = rng.randrange(n)
    cand = [w for w in V if w != u and frozenset((u, w)) not in Es]
    if not cand:
        return None
    heads = rng.sample(cand, rng.randint(1, min(heads_max, len(cand))))
    return V, E, [(u, w) for w in heads]


@pytest.mark.parametrize('seed', range(4))
def test_star_rule_matches_bruteforce_every_position(seed):
    rng = random.Random(seed)
    checked = 0
    for _ in range(60):
        inst = _random_instance(rng, rng.randint(3, 8), 4)
        if inst is None:
            continue
        V, E, arcs = inst
        u = arcs[0][0]
        for r in range(1, len(V) + 1):
            for K in itertools.combinations(V, r):
                for s in K:
                    assert gamma_outcome(V, E, arcs, s, K) == gamma_outcome_bruteforce(V, E, arcs, s, K)
                    checked += 1
    assert checked > 1000


def test_single_arc_is_theorem_10_11():
    # s = u: u essential in H + uw; s = w: w essential in H
    V, E = range(4), [(0, 1), (1, 2)]           # path 0-1-2 and isolated 3; arc 3 -> 0
    arcs = [(3, 0)]
    for s in V:
        assert gamma_outcome(V, E, arcs, s) == gamma_outcome_bruteforce(V, E, arcs, s)


def test_star_rule_rejects_two_tails():
    with pytest.raises(NotStarShaped):
        gamma_outcome(range(4), [], [(0, 2), (1, 3)], 0)


def test_dvg_solver_matches_bruteforce():
    rng = random.Random(3)
    digs = 0
    while digs < 60:
        n = rng.randint(3, 9)
        arcs = set()
        for a in range(n):
            for b in range(a + 1, n):
                r = rng.random()
                if r < 0.22:
                    arcs |= {(a, b), (b, a)}
                elif r < 0.34:
                    arcs.add((a, b) if rng.random() < 0.5 else (b, a))
        D = Digraph(n, arcs)
        try:
            solve_dvg(D, 0)
        except NotStarShaped:
            continue
        digs += 1
        win = D.outcome_function()
        for (W, v) in D.reachable_positions():
            assert solve_dvg(D, v, W) == win(W, v)


def test_barrier_examples_have_equal_key_local_data_but_opposite_outcomes():
    # Proposition 12.5 (disjoint arcs)
    K = _edges('0-8 0-9 1-8 1-9 2-6 2-7 2-8 3-7 3-8 4-6 4-8 5-6 5-7 5-9')
    A = [(1, 4), (2, 0)]
    assert gamma_outcome_bruteforce(range(10), K, A, 5) and not gamma_outcome_bruteforce(range(10), K, A, 3)
    assert key_local_data(range(10), K, [5, 1, 4, 2, 0]) == key_local_data(range(10), K, [3, 1, 4, 2, 0])
    # Proposition 12.6 (chained)
    K1, A1 = _edges('0-2 1-2 1-4 2-5 3-4 4-5'), [(0, 3), (3, 1)]
    K2 = _edges('0-6 0-7 0-8 0-9 1-6 1-7 1-9 2-6 2-8 2-9 3-6 3-8 3-9 4-6 4-7 4-9 5-6 5-7 5-9')
    A2 = [(1, 3), (3, 4)]
    assert gamma_outcome_bruteforce(range(6), K1, A1, 5) and not gamma_outcome_bruteforce(range(10), K2, A2, 0)
    assert key_local_data(range(6), K1, [5, 0, 3, 3, 1]) == key_local_data(range(10), K2, [0, 1, 3, 3, 4])
    # Proposition 12.9 (G6: s,h,h',a,b,l = 0..5; chain a->b->l versus l->b->a)
    G6 = _edges('0-1 1-2 1-3 1-4 2-3 2-4 2-5')
    assert gamma_outcome_bruteforce(range(6), G6, [(3, 4), (4, 5)], 0)
    assert not gamma_outcome_bruteforce(range(6), G6, [(5, 4), (4, 3)], 0)


def test_barrier_common_head_example():
    # Proposition 12.6 (common head): equal key-local data, opposite outcomes
    K1 = _edges('0-5 0-6 1-5 1-6 2-6 2-7 3-5 4-6 4-7')
    A1 = [(1, 4), (0, 4)]
    K2 = _edges('0-5 0-6 0-7 1-5 1-6 1-7 2-5 2-7 3-5 3-6 4-5 4-6 4-7')
    A2 = [(4, 0), (1, 0)]
    assert gamma_outcome_bruteforce(range(8), K1, A1, 3)
    assert not gamma_outcome_bruteforce(range(8), K2, A2, 2)
    assert key_local_data(range(8), K1, [3, 1, 4, 0, 4]) == key_local_data(range(8), K2, [2, 4, 0, 1, 0])
