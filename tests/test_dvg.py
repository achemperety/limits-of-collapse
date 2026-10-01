import random
from dvgcollapse.dvg import Digraph, random_digraph, order_at_zero
from dvgcollapse.potentials import check_outcomes_by_valuation


def test_directed_cycle_positions_and_resolvent():
    C5 = Digraph(5, [(i, (i + 1) % 5) for i in range(5)])
    # reachable positions: intervals of length 1..5 starting at their first vertex
    assert len(C5.reachable_positions()) == 5 * 5
    R = C5.resolvent_function()
    r = R(frozenset(range(5)), 0)          # a directed path on 5 vertices: mu(P4)/mu(P5)
    assert r.denom.degree() == 5 and order_at_zero(r) == -1   # P5 has odd order: P-position
    assert C5.wins_from(0) is False


def test_strong_components_and_scc_symmetry():
    D = Digraph(5, [(0, 1), (1, 0), (1, 2), (2, 3), (3, 2), (3, 4)])
    comps = sorted(sorted(C) for C in D.strong_components())
    assert comps == [[0, 1], [2, 3], [4]]
    assert D.is_scc_symmetric()
    T = Digraph(3, [(0, 1), (1, 0), (0, 2), (2, 0), (1, 2)])
    assert not T.is_scc_symmetric()
    assert T.one_way_arcs_in_components() == [(1, 2)]


def test_theorem_5_2_on_random_digraphs():
    rng = random.Random(11)
    for _ in range(25):
        D = random_digraph(rng.randint(3, 7), 0.4, rng)
        assert check_outcomes_by_valuation(D)


def test_gamma_and_back_arcs():
    T = Digraph(3, [(0, 1), (1, 2), (2, 0)])
    assert T.back_arcs([0, 1, 2]) == 1
    g = T.gamma([0, 1])
    # c(arc) = 1/(z^2 - 2) for the directed triangle (Theorem 11.3)
    assert str(g) == '1/(z**2 - 2)'
