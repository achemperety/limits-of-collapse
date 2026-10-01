"""Theorem 11.12 on five vertices: admissible self-energies of the 14 cores that pass the back-arc test."""
import pickle
from eta_variety import analyse
cores = pickle.load(open('cores5_backok.pkl', 'rb'))
for arcs in cores:
    kind, I = analyse(arcs)
    if kind != 'variety':
        print(kind.upper(), arcs, flush=True); continue
    print('arcs', arcs, 'dim', I.dimension(), flush=True)
    for P in I.minimal_associated_primes():
        print('    component dim', P.dimension(), ':', [str(g) for g in P.gens()], flush=True)
