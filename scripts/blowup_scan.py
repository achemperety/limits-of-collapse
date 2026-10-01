"""Conjecture 13.16(1): the continuant criterion (Theorem 13.9) on all size vectors with 3 <= m <= 8
and sizes at most 4, 4, 3, 3, 2, 2 (up to rotation); a set potential exists exactly for the
two-periodic ones.  Also: two-periodic families (s, t)^r for r <= 5, s, t <= 5."""
import itertools, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
from dvgcollapse.blowups import criterion, two_periodic

bad, cnt = [], 0
for m in range(3, 9):
    smax = {3: 4, 4: 4, 5: 3, 6: 3, 7: 2, 8: 2}[m]
    for sizes in itertools.product(range(1, smax + 1), repeat=m):
        if sizes != min(sizes[i:] + sizes[:i] for i in range(m)):
            continue
        cnt += 1
        if criterion(list(sizes)) != two_periodic(sizes):
            bad.append(sizes)
print('size vectors tested (up to rotation):', cnt, '; counterexamples to Conjecture 13.16(1):', bad)
ok = all(criterion([s, t] * r) for r in range(2, 6) for s in range(1, 6) for t in range(1, 6))
print('two-periodic (s,t)^r, r = 2..5, s,t <= 5: all have set potentials:', ok)
