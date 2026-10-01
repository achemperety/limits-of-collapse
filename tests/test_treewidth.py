from dvgcollapse.treewidth import ComponentStatePotential, ladder, resolvent_degree


def test_ladder_degrees_theorem_10_14():
    for m in range(2, 6):
        assert resolvent_degree(ladder(m), 0) == 2 ** (m + 1) - 2


def test_component_state_potential_agrees_with_game_tree():
    for m in range(2, 5):
        L = ladder(m)
        P = ComponentStatePotential(L)
        R = L.resolvent_function()
        assert all(P.resolvent(W, v) == R(W, v) for (W, v) in L.reachable_positions())
        assert P.number_of_states() == 8 * m


def test_resolvent_depends_only_on_canonical_subgame():
    import random
    from dvgcollapse.dvg import random_digraph
    from dvgcollapse.treewidth import canonical_subgame
    rng = random.Random(7)
    for _ in range(30):
        D = random_digraph(6, 0.35, rng)
        R = D.resolvent_function()
        seen = {}
        for (W, v) in D.reachable_positions():
            key = canonical_subgame(D, W, v)
            val = R(W, v)
            assert seen.setdefault(key, val) == val


def test_comb_has_exponentially_many_resolvents_at_the_collector():
    from dvgcollapse.treewidth import resolvents_at_collector
    for m in range(1, 5):
        assert len(resolvents_at_collector(m)) == 2 ** m


def test_closed_ladder_is_strongly_connected_and_not_symmetric():
    from dvgcollapse.treewidth import closed_ladder
    for m in range(2, 6):
        D = closed_ladder(m)
        assert len(D.strong_components()) == 1 and not D.is_scc_symmetric()
