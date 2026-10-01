"""Can the admissible self-energies of the bowtie and of the (1,2,2) triangle blow-up be produced by exits?
Library: root resolvents of all fresh games on at most 4 vertices; self-energies = sums of at most 3 of them."""
import itertools
import dvg
Kz = dvg.Kz; z = Kz(dvg.z)
lib = set()
for n in range(1, 5):
    for arcs in dvg.all_digraphs(n):
        R = dvg.make_resolvent(n, arcs)
        V = frozenset(range(n))
        reach = dvg.reach_matrix(n, arcs)
        for t in range(n):
            if all(reach[t][u] for u in range(n)):     # t reaches every vertex: a fresh game rooted at t
                lib.add(R(V, t))
lib = sorted(lib, key=str)
sums = {Kz(0)}
for r in (1, 2, 3):
    for combo in itertools.combinations_with_replacement(lib, r):
        sums.add(sum(combo, Kz(0)))
print('library of fresh root resolvents:', len(lib), '| exit self-energies (sums of <= 3):', len(sums))
bow = [e - 1 / ((z - e) - 1 / (z - e)) for e in sums]
print('bowtie, first component: realizable pairs found:', sum(1 for x in bow if x in sums))
blow = [a - 1 / (z - a) for a in sums]
print('(1,2,2) blow-up: realizable pairs found:', sum(1 for x in blow if x in sums))
book = [e + 1 / (z - e) for e in sums]
print('control, 2-page book: realizable pairs found:', sum(1 for x in book if x in sums))
