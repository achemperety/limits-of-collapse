"""The unified solver agrees with exhaustive search, whichever method it dispatches to."""
import random

from dvgcollapse.dvg import Digraph, random_digraph
from dvgcollapse.solver import main, outcome


def test_solver_matches_exhaustive_search():
    rng = random.Random(11)
    seen = set()
    for _ in range(150):
        n = rng.randint(2, 8)
        D = random_digraph(n, rng.uniform(0.2, 0.5), rng)
        win = D.outcome_function()
        for (W, v) in list(D.reachable_positions())[:60]:
            got, method = outcome(D, v, W)
            seen.add(method)
            assert got == win(W, v)
    assert {'common-tail', 'component-states'} <= seen


def test_exhaustive_fallback():
    D = Digraph(4, [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2), (1, 3)])
    got, method = outcome(D, 0, max_component=2)
    assert method == 'exhaustive' and got == D.wins_from(0)


def test_cli(capsys):
    assert main(['--n', '3', '--arcs', '0-1 1-2 2-0', '--start', '0']) == 0
    assert 'first player' in capsys.readouterr().out
