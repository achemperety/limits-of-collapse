# Verify the weighted two-point (Heilmann-Lieb / Godsil) identity with indeterminate vertex weights:
#   mu(S-x) mu(S-y) - mu(S) mu(S-x-y) = sum over U-paths P from x to y inside S of mu(S \ P)^2
import itertools, random
from sage.all__sagemath_singular import PolynomialRing, QQ

def mu(S, E, a, memo):
    S = frozenset(S)
    if S in memo: return memo[S]
    if not S: return a[0]**0
    v = min(S); rest = S - {v}
    val = a[v] * mu(rest, E, a, memo)
    for w in rest:
        if (min(v, w), max(v, w)) in E:
            val -= mu(rest - {w}, E, a, memo)
    memo[S] = val
    return val

def paths(x, y, S, E):
    out = []
    def rec(p):
        u = p[-1]
        if u == y: out.append(tuple(p)); return
        for w in S:
            if w not in p and (min(u, w), max(u, w)) in E:
                rec(p + [w])
    rec([x]); return out

random.seed(1)
bad = 0; tests = 0
for n in range(2, 8):
    R = PolynomialRing(QQ, ['a%d' % i for i in range(n)])
    a = R.gens()
    for trial in range(15):
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        E = set(p for p in pairs if random.random() < 0.5)
        memo = {}
        V = frozenset(range(n))
        for x, y in itertools.permutations(range(n), 2):
            lhs = mu(V - {x}, E, a, memo) * mu(V - {y}, E, a, memo) - mu(V, E, a, memo) * mu(V - {x, y}, E, a, memo)
            rhs = sum((mu(V - set(P), E, a, memo)**2 for P in paths(x, y, V, E)), R(0))
            tests += 1; bad += (lhs != rhs)
print("weighted two-point identity: tests", tests, "failures", bad)
