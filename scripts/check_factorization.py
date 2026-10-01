"""Symbolic check of the first-order factorization lemma.
Omega = prod_{t in T}(1+v_t) + eps*E(v), E multi-affine with symbolic coefficients lam[X].
Claim: d/d eps of Delta_xy(Omega) at eps=0 equals Q(v_R) * sum_S kappa(S) v^S,
Q = prod_{r in R}(1+v_r), kappa(S) = lam[S+x]+lam[S+y]-lam[S+x+y]-lam[S]."""
import itertools, sympy as sp
for n in range(2, 6):
    v = sp.symbols('v0:%d' % n)
    x, y = 0, 1
    R = list(range(2, n))
    subsets = [frozenset(c) for k in range(n + 1) for c in itertools.combinations(range(n), k)]
    lam = {X: sp.Symbol('l_' + ''.join(map(str, sorted(X))) if X else 'l_e') for X in subsets}
    eps = sp.Symbol('eps')
    Om = sp.prod([1 + v[t] for t in range(n)]) + eps * sum(lam[X] * sp.prod([v[t] for t in X]) for X in subsets)
    Om = sp.expand(Om)
    D = sp.diff(Om, v[x]) * sp.diff(Om, v[y]) - Om * sp.diff(Om, v[x], v[y])
    first = sp.expand(sp.diff(D, eps).subs(eps, 0))
    zero = sp.expand(D.subs(eps, 0))
    Q = sp.prod([1 + v[r] for r in R])
    F = 0
    for k in range(len(R) + 1):
        for c in itertools.combinations(R, k):
            S = frozenset(c)
            kap = lam[S | {x}] + lam[S | {y}] - lam[S | {x, y}] - lam[S]
            F += kap * sp.prod([v[r] for r in S])
    ok = sp.expand(first - Q * F) == 0
    print('n=%d  Delta(eps=0) == 0: %s   first-order == Q*F: %s' % (n, zero == 0, ok))
