"""Predictions of the structure theorems for set potentials, checked exactly over Q(z) on larger digraphs."""
from setpot_symbolic import set_potential_exact
import dvg
def sym(a, b): return [(a, b), (b, a)]
tests = []
# directed triangle 0->1->2->0
T = [(0, 1), (1, 2), (2, 0)]
tests.append(('triangle, a private sink at every vertex', 6, T + [(0, 3), (1, 4), (2, 5)], True))
tests.append(('triangle, sinks at two vertices only', 5, T + [(0, 3), (1, 4)], False))
tests.append(('triangle, pendant at 0, sinks at 1 and 2', 6, T + sym(0, 3) + [(1, 4), (2, 5)], True))
tests.append(('triangle, pendant at 0, sink at 1 only', 5, T + sym(0, 3) + [(1, 4)], False))
# hub with a pendant path of length 2 (resolvent z/(z^2-1)); others exit to a fresh directed 2-path t -> t'
tests.append(('triangle, pendant 2-path at 0, exits to directed 2-paths at 1 and 2', 9,
              T + sym(0, 3) + sym(3, 4) + [(1, 5), (5, 6), (2, 7), (7, 8)], True))
tests.append(('same, but vertex 2 exits to a sink instead', 8,
              T + sym(0, 3) + sym(3, 4) + [(1, 5), (5, 6), (2, 7)], False))
C4 = [(0, 1), (1, 2), (2, 3), (3, 0)]
tests.append(('4-cycle, sinks at opposite vertices', 6, C4 + [(0, 4), (2, 5)], True))
tests.append(('4-cycle, sinks at adjacent vertices', 6, C4 + [(0, 4), (1, 5)], False))
tests.append(('4-cycle, pendant at 0 and sink at 2', 6, C4 + sym(0, 4) + [(2, 5)], True))
# books: spine x=0 -> y=1, pages 2..r+1 with 1 -> p -> 0
def book(r): return [(0, 1)] + [(1, p) for p in range(2, r + 2)] + [(p, 0) for p in range(2, r + 2)]
tests.append(('book with 3 pages, two pendants at the spine tail', 7, book(3) + sym(0, 5) + sym(0, 6), True))
tests.append(('book with 3 pages, one pendant at the spine tail', 6, book(3) + sym(0, 5), False))
tests.append(('book with 4 pages, three sinks at the spine tail', 9, book(4) + [(0, 6), (0, 7), (0, 8)], True))
tests.append(('book with 2 pages, pendant at tail, a sink at every other vertex and 2 at tail?', 9,
              book(2) + sym(0, 4) + [(1, 5), (2, 6), (3, 7), (0, 8)], False))
# 5-cycle with 2-periodic self-energies is impossible (odd): sinks at 0 and 2 only
C5 = [(i, (i + 1) % 5) for i in range(5)]
tests.append(('5-cycle, sinks at 0 and 2', 7, C5 + [(0, 5), (2, 6)], False))
tests.append(('5-cycle, a sink at every vertex', 10, C5 + [(i, 5 + i) for i in range(5)], True))
# 6-cycle with 2-periodic exits
C6 = [(i, (i + 1) % 6) for i in range(6)]
tests.append(('6-cycle, sinks at 0, 2, 4', 9, C6 + [(0, 6), (2, 7), (4, 8)], True))
tests.append(('6-cycle, sinks at 0, 2, 3', 9, C6 + [(0, 6), (2, 7), (3, 8)], False))
ok = 0
for name, n, arcs, pred in tests:
    got = set_potential_exact(n, arcs)
    ok += (got == pred)
    print('%-80s predicted %-5s exact %-5s %s' % (name, pred, got, '' if got == pred else '<-- MISMATCH'))
print(ok, 'of', len(tests), 'predictions confirmed')
