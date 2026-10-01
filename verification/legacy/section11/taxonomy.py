"""Frontier III: taxonomy of non-SCC-symmetric digraphs with a set potential (n <= 5), and a check of the
component reduction: D has a set potential iff every non-symmetric strong component C, played with its exit
self-energies eta_x = sum of fresh root resolvents of the exits of x, has one."""
import pickle, itertools
from functools import lru_cache
from collections import Counter, defaultdict
import dvg
from setpot_symbolic import set_potential_exact

Kz, z = dvg.Kz, dvg.Kz(dvg.z)


def sccs(n, arcs):
    reach = dvg.reach_matrix(n, arcs)
    comps, seen = [], set()
    for v in range(n):
        if v in seen:
            continue
        C = frozenset(u for u in range(n) if reach[v][u] and reach[u][v])
        seen |= C
        comps.append(C)
    return comps, reach


def fresh_resolvent(n, arcs):
    """rho(t) = R(Reach(t), t): root resolvent of the fresh game at t (every vertex unvisited)."""
    R = dvg.make_resolvent(n, arcs)
    V = frozenset(range(n))
    return lambda t: R(V, t)


def component_game_setpot(C, arcs, eta):
    """Set potential of the game on D[C] with self-energies eta (exact over Q(z))."""
    out = {v: [w for (a, w) in arcs if a == v and w in C] for v in C}

    @lru_cache(maxsize=None)
    def R(W, v):
        s = eta[v]
        for w in out[v]:
            if w in W and w != v:
                s += R(W - {v}, w)
        return 1 / (z - s)
    # reachable positions inside C
    seen = set()
    stack = [(frozenset(C), s) for s in C]
    while stack:
        W, v = stack.pop()
        if (W, v) in seen:
            continue
        seen.add((W, v))
        for w in out[v]:
            if w in W and w != v:
                stack.append((W - {v}, w))
    M = {frozenset(C): Kz(1)}
    for (W, v) in sorted(seen, key=lambda p: -len(p[0])):
        if W not in M:
            continue
        m2 = M[W] * R(W, v)
        W2 = W - {v}
        if W2 in M:
            if M[W2] != m2:
                return False
        else:
            M[W2] = m2
    return True


def is_directed_cycle(C, arcs):
    A = [(a, b) for (a, b) in arcs if a in C and b in C]
    if len(A) != len(C):
        return False
    outd = Counter(a for a, b in A); ind = Counter(b for a, b in A)
    return all(outd[v] == 1 and ind[v] == 1 for v in C) and len(C) >= 3


def cycle_order(C, arcs):
    nxt = {a: b for (a, b) in arcs if a in C and b in C}
    start = min(C); order = [start]
    while nxt[order[-1]] != start:
        order.append(nxt[order[-1]])
    return order


def describe(n, arcs):
    comps, reach = sccs(n, arcs)
    A = set(arcs)
    rho = fresh_resolvent(n, arcs)
    info = []
    for C in comps:
        sym = all((b, a) in A for (a, b) in A if a in C and b in C)
        if sym:
            continue
        eta = {x: sum((rho(t) for t in range(n) if (x, t) in A and t not in C), Kz(0)) for x in C}
        cyc = is_directed_cycle(C, arcs)
        entry = sorted(x for x in C if any((y, x) in A for y in range(n) if y not in C))
        d = {'size': len(C), 'cycle': cyc, 'eta': eta, 'entry': entry,
             'canon': dvg.canon(len(C), [(sorted(C).index(a), sorted(C).index(b)) for (a, b) in arcs
                                         if a in C and b in C])}
        if cyc:
            order = cycle_order(C, arcs)
            d['eta_seq'] = [eta[x] for x in order]
            d['order'] = order
        d['comp_setpot'] = component_game_setpot(C, arcs, eta)
        info.append(d)
    return comps, info


if __name__ == '__main__':
    import dvg as _d
    census = {}
    allD = {3: list(_d.all_digraphs(3)), 4: list(_d.all_digraphs(4)), 5: pickle.load(open('digraphs5.pkl', 'rb'))}
    mismatch = 0
    rows = []
    for n in (3, 4, 5):
        for arcs in allD[n]:
            if _d.scc_symmetric(n, arcs):
                continue
            glob = set_potential_exact(n, arcs)
            comps, info = describe(n, arcs)
            loc = all(d['comp_setpot'] for d in info)
            if glob != loc:
                mismatch += 1
                print('MISMATCH', n, arcs, glob, loc)
            if glob:
                rows.append((n, arcs, comps, info))
    print('component reduction mismatches:', mismatch)
    print('non-SCC-symmetric digraphs with a set potential:', Counter(r[0] for r in rows))
    # taxonomy by the multiset of non-symmetric components
    tax = defaultdict(list)
    for n, arcs, comps, info in rows:
        key = tuple(sorted((d['size'], d['cycle'], d['canon']) for d in info))
        tax[(n, key)].append(arcs)
    for (n, key), lst in sorted(tax.items()):
        print(n, [(s, 'cycle' if c else 'NONCYCLE', can) for s, c, can in key], len(lst))
    pickle.dump(rows, open('setpot_rows.pkl', 'wb'))
