"""Census of non-SCC-symmetric digraphs with a set potential, from the back-arc survivors (exact test over Q(z)),
and the pruned cores of their non-symmetric strong components."""
import sys, time, pickle
from collections import Counter
import dvg
from setpot_symbolic import set_potential_exact
from describe_cores import parse, symmetric_branches, blowup


def components(n, arcs):
    reach = dvg.reach_matrix(n, arcs)
    comps, seen = [], set()
    for v in range(n):
        if v in seen:
            continue
        C = frozenset(u for u in range(n) if reach[v][u] and reach[u][v])
        seen |= C; comps.append(C)
    return comps


def prune(C, arcs):
    """Remove symmetric branches repeatedly; return the pruned vertex set and the induced arcs (relabelled)."""
    C = set(C)
    while True:
        verts = sorted(C)
        idx = {v: i for i, v in enumerate(verts)}
        sub = [(idx[a], idx[b]) for (a, b) in arcs if a in C and b in C]
        br = symmetric_branches(len(verts), sub)
        if not br:
            return len(verts), sub
        h, comp = br[0]
        C -= {verts[i] for i in comp}


if __name__ == '__main__':
    fname = sys.argv[1]
    rows = [parse(L) for L in open(fname) if L.strip()]
    t0 = time.time(); found = []
    for i, (n, arcs) in enumerate(rows):
        if set_potential_exact(n, arcs):
            found.append((n, arcs))
        if i % 1000 == 0:
            print(i, len(found), round(time.time() - t0), flush=True)
    kinds = Counter()
    for n, arcs in found:
        A = set(arcs)
        for C in components(n, arcs):
            if all((b, a) in A for (a, b) in A if a in C and b in C):
                continue
            k, sub = prune(C, arcs)
            b = blowup(k, sub)
            kinds[(len(C), k, str(min(b[i:] + b[:i] for i in range(len(b)))) if b else 'not a blow-up: ' + str(sub))] += 1
    print('with a set potential:', len(found), 'of', len(rows), round(time.time() - t0), 's')
    for key, v in sorted(kinds.items()):
        print('   component size %d, pruned core size %d, %s: %d' % (key[0], key[1], key[2], v))
    pickle.dump(found, open(fname + '.setpot.pkl', 'wb'))
