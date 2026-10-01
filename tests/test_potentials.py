import random
from dvgcollapse.dvg import Digraph
from dvgcollapse.potentials import scc_symmetric_potential, verify_potential, set_potential, valuation_potential


def test_theorem_6_6_potential_on_random_scc_symmetric_digraphs():
    rng = random.Random(5)
    done = 0
    while done < 25:
        n = rng.randint(3, 8)
        arcs = set()
        for a in range(n):
            for b in range(a + 1, n):
                r = rng.random()
                if r < 0.25:
                    arcs |= {(a, b), (b, a)}
                elif r < 0.45:
                    arcs.add((a, b))
        D = Digraph(n, arcs)
        if not D.is_scc_symmetric():
            continue
        U, a = scc_symmetric_potential(D)
        assert verify_potential(D, U, a)
        done += 1


def test_hierarchy_examples_of_theorem_10_13():
    C4 = Digraph(4, [(0, 1), (1, 2), (2, 3), (3, 0)])
    assert set_potential(C4) is not None                         # directed cycle: set potential
    D = Digraph(4, [(0, 1), (1, 0), (0, 2), (2, 1), (2, 3)])
    assert set_potential(D) is None and valuation_potential(D) is not None
    T = Digraph(3, [(0, 1), (1, 0), (0, 2), (2, 0), (1, 2)])
    assert set_potential(T) is None and valuation_potential(T) is None
