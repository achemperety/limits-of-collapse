"""The complex counterexample to Theorem 6.10 over the algebraic closure of Q(z) (Section 11).

For m >= 3 and r >= 0 let D(m, r) be the digraph with

  * a directed m-cycle c_0 -> c_1 -> ... -> c_{m-1} -> c_0,
  * r private sinks at every cycle vertex (arcs c_i -> s_{i,1..r}),
  * m + 1 isolated vertices f_0, ..., f_m (the environment F).

D(m, 0) is D_m of Theorem 11.6.  The potential (U, a) is

  * U = all edges between the cycle and F (a complete bipartite graph K_{m, m+1}),
  * a = 0 on the cycle, a = z on every sink, and on F the m + 1 roots of
        Y(y) = y^(m+1) - (m! z / mu(P_m)(x)) * sum_{l=0}^{m} mu(P_l)(x) y^l / l!,   x = z - r/z,
    where mu(P_l) is the matching polynomial of the path with l vertices.

Every value mu_a(U[S]) that the potential identities involve is a symmetric polynomial in the
environment weights, so each identity can be verified exactly by rewriting it in the elementary
symmetric functions e_1, ..., e_{m+1}, which Vieta's formulas fix from Y.  `verify_symbolic`
does this for every reachable position of D(m, r), with resolvents computed independently by the
game recursion of `Digraph.resolvent_function`.  `verify_numeric` evaluates the same identities
with the numerically computed roots of Y at a (complex) value of z.
"""
from __future__ import annotations

import itertools
from typing import Dict, List, Tuple

import sympy as sp

from .dvg import Digraph

z = sp.Symbol('z')
y = sp.Symbol('y')


def mu_path(l: int, x=z):
    """Matching polynomial mu(P_l)(x) of the path with l vertices."""
    a, b = sp.Integer(1), x
    if l == 0:
        return a
    for _ in range(l - 1):
        a, b = b, sp.expand(x * b - a)
    return b


def counterexample(m: int, r: int = 0) -> Tuple[Digraph, List[Tuple[int, int]], Dict[str, List[int]]]:
    """Return (D, U, parts) for D(m, r); parts has keys 'cycle', 'env', 'sinks'."""
    if m < 3:
        raise ValueError('m >= 3 is required (a 2-cycle is a symmetric pair)')
    cycle = list(range(m))
    sinks = list(range(m, m + m * r))
    env = list(range(m + m * r, m + m * r + m + 1))
    arcs = [(i, (i + 1) % m) for i in range(m)]
    for i in range(m):
        for j in range(r):
            arcs.append((i, m + i * r + j))
    D = Digraph(len(cycle) + len(sinks) + len(env), arcs)
    U = [(c, f) for c in cycle for f in env]
    return D, U, {'cycle': cycle, 'env': env, 'sinks': sinks}


def prescribed_elementary(m: int, r: int = 0) -> Dict[int, sp.Expr]:
    """e_1, ..., e_{m+1} of the environment weights: e_{m+1-l} = (-1)^l c mu(P_l)(x) / l!."""
    x = z - sp.Rational(r) / z if r else z
    c = (-1) ** m * sp.factorial(m) * z / mu_path(m, x)
    return {m + 1 - l: sp.simplify((-1) ** l * c * mu_path(l, x) / sp.factorial(l)) for l in range(m + 1)}


def Y_polynomial(m: int, r: int = 0) -> sp.Expr:
    """The monic polynomial whose roots are the environment weights (in y, coefficients in Q(z))."""
    e = prescribed_elementary(m, r)
    return sp.together(y ** (m + 1) + sum((-1) ** j * e[j] * y ** (m + 1 - j) for j in range(1, m + 2)))


def Y_cleared(m: int, r: int = 0) -> sp.Expr:
    """Y_m with denominators cleared and content removed (a polynomial in Z[z, y])."""
    num = sp.numer(sp.together(Y_polynomial(m, r)))
    P = sp.Poly(sp.expand(num), y, z)
    content = sp.gcd_list([sp.Integer(c) for c in P.coeffs()])
    expr = sp.expand(num / content)
    g = sp.gcd_list(sp.Poly(expr, y).all_coeffs())
    return sp.expand(sp.cancel(expr / g))


def _weights_symbolic(m: int, r: int, parts) -> Tuple[Dict[int, sp.Expr], Tuple[sp.Symbol, ...]]:
    beta = sp.symbols('b0:%d' % (m + 1))
    a: Dict[int, sp.Expr] = {c: sp.Integer(0) for c in parts['cycle']}
    a.update({s: z for s in parts['sinks']})
    a.update({f: beta[i] for i, f in enumerate(parts['env'])})
    return a, beta


def _mu_factory(U, a):
    adj: Dict[int, set] = {}
    for (i, j) in U:
        adj.setdefault(i, set()).add(j)
        adj.setdefault(j, set()).add(i)
    memo: Dict[frozenset, sp.Expr] = {}

    def mu(S):
        S = frozenset(S)
        if S in memo:
            return memo[S]
        if not S:
            return sp.Integer(1)
        v = min(S)
        rest = S - {v}
        val = a[v] * mu(rest)
        for w in adj.get(v, ()):
            if w in rest:
                val -= mu(rest - {w})
        memo[S] = sp.expand(val)
        return memo[S]

    return mu


def verify_symbolic(m: int, r: int = 0) -> Dict[str, object]:
    """Exact verification of every potential identity of D(m, r) over Q(z)(e_1..e_{m+1}).

    Returns a dictionary with the number of reachable positions checked and a Boolean `ok`."""
    D, U, parts = counterexample(m, r)
    a, beta = _weights_symbolic(m, r, parts)
    mu = _mu_factory(U, a)
    e = prescribed_elementary(m, r)
    R = D.resolvent_function()

    cache: Dict[frozenset, sp.Expr] = {}

    def value(S):
        S = frozenset(S)
        if S not in cache:
            poly = mu(S)
            if poly.free_symbols & set(beta):
                sym, rem, defs = sp.polys.polyfuncs.symmetrize(poly, *beta, formal=True)
                if rem != 0:
                    raise AssertionError('mu_a(U[S]) is not symmetric in the environment weights')
                sub = {s: e[int(str(s)[1:])] for s, _ in defs}
                cache[S] = sp.simplify(sym.subs(sub))
            else:
                cache[S] = sp.simplify(poly)
        return cache[S]

    positions = D.reachable_positions()
    ok = True
    for (W, v) in positions:
        den = value(W)
        num = value(W - {v})
        Rexpr = R(W, v).as_expr()
        if sp.simplify(den) == 0 or sp.simplify(num - Rexpr * den) != 0:
            ok = False
            break
    return {'m': m, 'r': r, 'positions': len(positions), 'ok': ok}


def verify_numeric(m: int, zval, r: int = 0, dps: int = 50) -> Dict[str, object]:
    """Numerical verification at a real or complex z: roots of Y, then every identity."""
    import mpmath as mp
    mp.mp.dps = dps
    zv = mp.mpmathify(zval)
    D, U, parts = counterexample(m, r)
    Yexpr = sp.Poly(Y_cleared(m, r), y)
    coeffs = [mp.mpmathify(sp.lambdify(z, c, 'mpmath')(zv)) for c in Yexpr.all_coeffs()]
    roots = mp.polyroots(coeffs, maxsteps=500, extraprec=4 * dps)
    a = {c: mp.mpf(0) for c in parts['cycle']}
    a.update({s: zv for s in parts['sinks']})
    a.update({f: roots[i] for i, f in enumerate(parts['env'])})
    adj: Dict[int, set] = {}
    for (i, j) in U:
        adj.setdefault(i, set()).add(j)
        adj.setdefault(j, set()).add(i)
    memo: Dict[frozenset, object] = {}

    def mu(S):
        S = frozenset(S)
        if S in memo:
            return memo[S]
        if not S:
            return mp.mpf(1)
        v = min(S)
        rest = S - {v}
        val = a[v] * mu(rest)
        for w in adj.get(v, ()):
            if w in rest:
                val -= mu(rest - {w})
        memo[S] = val
        return val

    out = D.out
    rmemo: Dict[Tuple[frozenset, int], object] = {}

    def Rnum(W, v):
        key = (W, v)
        if key not in rmemo:
            W2 = W - {v}
            s = mp.mpf(0)
            for w in out[v]:
                if w in W2:
                    s += Rnum(W2, w)
            rmemo[key] = 1 / (zv - s)
        return rmemo[key]

    worst = mp.mpf(0)
    smallest = None
    for (W, v) in D.reachable_positions():
        den = mu(W)
        smallest = abs(den) if smallest is None else min(smallest, abs(den))
        res = abs(mu(W - {v}) / den - Rnum(W, v))
        worst = max(worst, res)
    nonreal = sum(1 for b in roots if abs(mp.im(b)) > mp.mpf(10) ** (-dps // 2))
    return {'m': m, 'r': r, 'z': zval, 'environment_weights': roots, 'max_residual': worst,
            'min_abs_denominator': smallest, 'nonreal_weights': nonreal}


def is_scc_symmetric(m: int, r: int = 0) -> bool:
    D, _, _ = counterexample(m, r)
    return D.is_scc_symmetric()
