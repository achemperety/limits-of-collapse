"""Smallest frozen environment for the directed 4-cycle (and 5-cycle with one frozen vertex):
does any graph U on C_m + F (|F| = f) admit a frozen potential over the algebraic closure of Q(z)?"""
import sys, time, itertools
from env_test import has_env_potential
m, f = int(sys.argv[1]), int(sys.argv[2])
arcs = [(i, (i + 1) % m) for i in range(m)]
N = m + f
pairs = [(a, b) for a in range(N) for b in range(a + 1, N)]
def images(U):
    for r in range(m):
        for perm in itertools.permutations(range(m, N)):
            mp_ = {i: (i + r) % m for i in range(m)}
            mp_.update({m + j: perm[j] for j in range(f)})
            yield tuple(sorted(tuple(sorted((mp_[a], mp_[b]))) for a, b in U))
seen = set(); found = []; t0 = time.time(); tested = 0
for mask in range(1 << len(pairs)):
    U = [pairs[i] for i in range(len(pairs)) if mask >> i & 1]
    c = min(images(U))
    if c in seen: continue
    seen.add(c)
    tested += 1
    if has_env_potential(m, arcs, f, frozenset(U)):
        found.append(U); print('FOUND', U, flush=True)
print('m', m, 'f', f, 'orbits tested', tested, 'with a frozen potential', len(found), round(time.time() - t0), 's')
