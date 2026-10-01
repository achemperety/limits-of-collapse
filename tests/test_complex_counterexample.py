"""Section 11: the counterexample to Theorem 6.10 over the algebraic closure of Q(z)."""
import sympy as sp

from dvgcollapse import complex_counterexample as cc


def test_counterexample_is_not_scc_symmetric():
    for m in (3, 4, 5):
        assert not cc.is_scc_symmetric(m)


def test_all_identities_exact_small_cases():
    for m, r in [(3, 0), (4, 0), (3, 1)]:
        res = cc.verify_symbolic(m, r)
        assert res['ok'], res


def test_Y3_explicit_form():
    z, y = cc.z, cc.y
    expected = (z**2 - 2) * y**4 - z * (z**2 - 2) * y**3 - 3 * (z**2 - 1) * y**2 - 6 * z * y - 6
    got = cc.Y_cleared(3)
    assert sp.expand(got - expected) == 0 or sp.expand(got + expected) == 0


def test_numeric_residuals_and_nonreal_weights_at_large_z():
    for m in (3, 4, 5):
        res = cc.verify_numeric(m, 40, dps=40)
        assert res['max_residual'] < 1e-30
        # Proposition 11.9(2): m-1 non-real weights for odd m, m for even m
        assert res['nonreal_weights'] == (m - 1 if m % 2 else m)
