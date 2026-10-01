"""Tests for Section 13: the continuant criterion for blow-ups and the twin identity."""
import itertools

from dvgcollapse import blowups as bu
from dvgcollapse.dvg import Z
from dvgcollapse.potentials import set_potential


def _rotation_classes(max_vertices):
    for m in range(3, max_vertices + 1):
        for sizes in itertools.product(range(1, max_vertices - m + 2), repeat=m):
            if sum(sizes) > max_vertices:
                continue
            if sizes != min(sizes[i:] + sizes[:i] for i in range(m)):
                continue
            yield sizes


def test_gamma_formula_on_all_paths():
    for sizes in [(2, 1, 1), (1, 2, 2), (2, 1, 2, 1), (2, 1, 1, 1), (3, 1, 1), (2, 2, 2), (1, 3, 1, 3)]:
        D = bu.blowup_digraph(sizes)
        R = D.resolvent_function()
        cls = [i for i, s in enumerate(sizes) for _ in range(s)]
        V = frozenset(range(D.n))

        def dfs(path, g):
            assert g == bu.gamma_formula(sizes, cls[path[0]], len(path))
            for w in D.out[path[-1]]:
                if w not in path:
                    dfs(path + [w], g * R(V - frozenset(path), w))

        for s0 in range(D.n):
            dfs([s0], R(V, s0))


def test_criterion_matches_direct_test_up_to_nine_vertices():
    count = 0
    for sizes in _rotation_classes(9):
        D = bu.blowup_digraph(sizes)
        direct = set_potential(D) is not None
        assert direct == bu.criterion(sizes), sizes
        assert direct == bu.two_periodic(sizes), sizes      # Conjecture 13.16(1) and Corollary 13.11
        count += 1
    assert count > 60


def test_twin_identity_symbolic():
    for m in (4, 6, 8, 10):
        for t in range(1, 4):
            for s in range(t + 1, t + 4):
                assert bu.twin_identity_holds(s, t, m)
    assert bu.twin_identity_holds(3, 1, 4, e=1 / Z)


def test_two_periodic_have_set_potentials():
    for r in range(2, 5):
        for s in range(1, 5):
            for t in range(1, 5):
                assert bu.criterion([s, t] * r)
