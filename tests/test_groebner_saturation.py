import pytest
from dvgcollapse.dvg import Digraph, all_digraphs
from dvgcollapse import groebner_saturation as gs

BACKENDS = ['sympy'] + (['singular'] if gs.have_singular() else [])


@pytest.mark.parametrize('backend', BACKENDS)
def test_three_vertex_census(backend):
    stats = {}
    for D in all_digraphs(3):
        key = (D.is_scc_symmetric(), bool(gs.search(D, backend=backend)))
        stats[key] = stats.get(key, 0) + 1
    assert stats == {(True, True): 13, (False, False): 3}


@pytest.mark.parametrize('backend', BACKENDS)
def test_frozen_triangle_complex_rescue(backend):
    C3 = Digraph(3, [(0, 1), (1, 2), (2, 0)])
    U = [(i, j) for i in range(3) for j in (3, 4)]          # the six edges between C and F
    assert gs.decide(C3, U, nf=2, backend=backend)           # a potential over the algebraic closure
    assert not gs.decide(C3, [(0, 3), (1, 3), (2, 3)], nf=1, backend=backend)


@pytest.mark.skipif(not gs.have_singular(), reason='Singular backend not installed')
def test_one_frozen_vertex_never_rescues_the_triangle():
    C3 = Digraph(3, [(0, 1), (1, 2), (2, 0)])
    assert gs.search(C3, nf=1, backend='singular') == []
