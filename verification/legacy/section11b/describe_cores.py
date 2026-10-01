"""Describe pruned cores: is the digraph a blow-up of a directed cycle (twin classes S_0 -> S_1 -> ... -> S_0)?"""
import sys
from collections import Counter


def parse(line):
    a = list(map(int, line.split())); nv, ne = a[0], a[1]
    return nv, [(a[2 + 2 * i], a[3 + 2 * i]) for i in range(ne)]


def symmetric_branches(nv, arcs):
    A = set(arcs); und = {v: set() for v in range(nv)}
    for a, b in arcs:
        und[a].add(b); und[b].add(a)
    out = []
    for h in range(nv):
        seen = {h}
        for s in range(nv):
            if s in seen:
                continue
            comp = set(); stack = [s]
            while stack:
                x = stack.pop()
                if x in seen:
                    continue
                seen.add(x); comp.add(x); stack.extend(und[x] - seen)
            S = comp | {h}
            if all((b, a) in A for (a, b) in A if a in S and b in S):
                out.append((h, sorted(comp)))
    return out


def blowup(nv, arcs):
    """Return the cyclic sequence of class sizes if the digraph is a blow-up of a directed cycle, else None."""
    A = set(arcs)
    inn = {v: frozenset(u for u in range(nv) if (u, v) in A) for v in range(nv)}
    out = {v: frozenset(u for u in range(nv) if (v, u) in A) for v in range(nv)}
    cls = {}
    for v in range(nv):
        cls.setdefault((inn[v], out[v]), []).append(v)
    classes = list(cls.values())
    if any((a, b) in A for c in classes for a in c for b in c):
        return None
    idx = {v: i for i, c in enumerate(classes) for v in c}
    succ = {}
    for i, c in enumerate(classes):
        targets = {idx[u] for u in out[c[0]]}
        if len(targets) != 1:
            return None
        t = targets.pop()
        if set(out[c[0]]) != set(classes[t]):
            return None
        succ[i] = t
    order = [0]
    while succ[order[-1]] != 0:
        order.append(succ[order[-1]])
        if len(order) > len(classes):
            return None
    if len(order) != len(classes):
        return None
    return tuple(len(classes[j]) for j in order)


if __name__ == '__main__':
    fname = sys.argv[1]
    cores = [parse(L) for L in open(fname) if L.strip()]
    pruned = [(nv, a) for nv, a in cores if not symmetric_branches(nv, a)]
    kinds = Counter()
    others = []
    for nv, a in pruned:
        b = blowup(nv, a)
        if b:
            rots = [b[i:] + b[:i] for i in range(len(b))]
            kinds[('blow-up of C%d' % len(b), min(rots))] += 1
        else:
            others.append((nv, a))
    print('back-arc survivors:', len(cores), '| pruned cores:', len(pruned))
    for k, v in sorted(kinds.items()):
        print('  ', k, v)
    print('  other pruned cores:', len(others))
    for nv, a in others:
        print('     ', a)
