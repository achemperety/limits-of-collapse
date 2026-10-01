"""Strongly connected non-symmetric 5-vertex cores that pass the back-arc test (Theorem 10.6):
all directed paths with the same vertex set have the same number of back arcs."""
import pickle
import dvg
digs = pickle.load(open('digraphs5.pkl', 'rb'))
n = 5
def strongly_connected(arcs):
    r = dvg.reach_matrix(n, arcs); return all(r[a][b] for a in range(n) for b in range(n))
def back_ok(arcs):
    A = set(arcs); out = dvg.out_nbrs(n, arcs)
    seen = {}
    def ext(path):
        X = frozenset(path)
        b = sum(1 for i in range(len(path)) for j in range(i) if (path[i], path[j]) in A)
        if seen.setdefault(X, b) != b: return False
        for w in out[path[-1]]:
            if w not in X and not ext(path + [w]): return False
        return True
    return all(ext([s]) for s in range(n))
cores = [a for a in digs if strongly_connected(a) and not dvg.scc_symmetric(n, a)]
ok = [a for a in cores if back_ok(a)]
print('strongly connected non-symmetric:', len(cores), '| back-arc invariant:', len(ok))
pickle.dump(ok, open('cores5_backok.pkl', 'wb'))
