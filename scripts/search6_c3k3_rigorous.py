"""Proposition 11.10: saturated Groebner test over Q(z) for all 1,928 graphs U (about 15 minutes)."""
import os
import sys, time, pickle
sys.argv = ['x']
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'search6_c3k3.py')).read().split('t0 = time.time(); surv = []')[0]
exec(src)
cand = [U for U in reps if triangle_connected(U)]
t0 = time.time(); found = []
for idx, U in enumerate(cand):
    if gs.decide(D, U, backend='singular'): found.append(U)
    if idx % 200 == 0: print(idx, len(found), round(time.time() - t0), 's', flush=True)
print('RIGOROUS RESULT over Q(z): C3 + 3 isolated vertices, %d graphs U (Theorem 10.7 filter applied): potentials found: %d %s (%d s)'
      % (len(cand), len(found), found, round(time.time() - t0)), flush=True)
