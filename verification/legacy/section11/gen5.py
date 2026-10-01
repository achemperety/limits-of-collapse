import numpy as np, itertools, pickle, time
n = 5
pairs = [(a, b) for a in range(n) for b in range(n) if a != b]
idx = {p: i for i, p in enumerate(pairs)}
N = 1 << len(pairs)
masks = np.arange(N, dtype=np.int64)
best = masks.copy()
t0 = time.time()
for perm in itertools.permutations(range(n)):
    newm = np.zeros(N, dtype=np.int64)
    for i, (a, b) in enumerate(pairs):
        j = idx[(perm[a], perm[b])]
        newm |= ((masks >> i) & 1) << j
    np.minimum(best, newm, out=best)
reps = np.unique(best)
print(len(reps), time.time() - t0)
digs = [[pairs[i] for i in range(len(pairs)) if (int(m) >> i) & 1] for m in reps]
pickle.dump(digs, open('digraphs5.pkl', 'wb'))
