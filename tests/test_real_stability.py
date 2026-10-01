from fractions import Fraction
from dvgcollapse.dvg import Digraph, all_digraphs
from dvgcollapse.real_stability import (first_order_factorization_holds, obstruction_certificate,
                                        rayleigh_difference, triangle_band_value)


def test_factorization_lemma_symbolic():
    for k in range(2, 6):
        assert first_order_factorization_holds(k)


def test_certificates_on_all_small_non_scc_symmetric_digraphs():
    count = 0
    for n in (3, 4):
        for D in all_digraphs(n):
            c = obstruction_certificate(D)
            if c is None:
                assert D.is_scc_symmetric()
                continue
            assert (c[3], c[4]) == (0, 1)
            count += 1
    assert count == 3 + 114


def test_exact_certificate_from_resolvents():
    D = Digraph(5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0), (2, 0), (3, 1)])
    assert obstruction_certificate(D, exact=True)[3:] == (0, 1)


def test_triangle_band_is_sharp():
    assert triangle_band_value(Fraction(19, 10)) > 0
    assert triangle_band_value(2) == 0
    assert triangle_band_value(Fraction(21, 10)) < 0 and triangle_band_value(50) < 0


def test_rayleigh_difference_of_a_product_vanishes():
    coeffs = {frozenset(s): 1.0 for s in [(), (0,), (1,), (2,), (0, 1), (0, 2), (1, 2), (0, 1, 2)]}
    assert abs(rayleigh_difference(coeffs, 0, 1, {2: -3.7})) < 1e-12
