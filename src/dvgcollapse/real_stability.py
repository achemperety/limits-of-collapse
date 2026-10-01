"""Real stability and the proof that real potentials force SCC-symmetry (Theorem 9.5 of the paper).

For a potential with real weights, g(h) = mu_{a+h}(U)/mu_a(U) = sum_T c(T) h^T is real stable
(Heilmann-Lieb), and c(T) = Gamma(P) for the vertex set T of any directed path P.  With
gamma_T = z^{|T|} c(T), Lemma 9.7 gives gamma_T = 1 + lambda_T z^{-2} + O(z^{-3}), where
    lambda_T = sum_{t in T} (d^+(t) + e_t) - back(P)                         (e_t: exits' weight)
for path sets, and the Tameness Lemma forces the same form on every subset of a directed cycle.
For a one-way arc x -> y on a directed cycle y = q_0 -> ... -> q_k = x -> y, with R = {q_1..q_{k-1}},
the first-order Rayleigh difference of the pair (x, y) is z^{-2} Q F with Q = prod_{r in R}(1 + v_r)
and F = sum_{S subset R} kappa(S) v^S, where
    kappa(S) = lambda_{S+x} + lambda_{S+y} - lambda_{S+x+y} - lambda_S.
The certificate below computes kappa(empty) = [y -> x] = 0 and kappa(R) = [x -> y] = 1 exactly
from path data; F = cQ is then impossible, which contradicts Brandén's criterion.

Functions
  rayleigh_difference(coeffs, i, j, point)   Delta_ij of a multi-affine polynomial at a real point
  lambda_of_path(D, path, exits)             the z^{-2} coefficient of gamma for a path set
  lambda_exact(D, path, eta)                 the same, read off the exact Gamma(P) in Q(z)
  obstruction_certificate(D)                 (x, y, cycle path, kappa(empty), kappa(R)) or None
  first_order_factorization_holds(k)         symbolic check of the factorization lemma, k variables
  triangle_band_value(z)                     Delta at v = -1/2 for the directed triangle (exact)
"""
from __future__ import annotations

import itertools
from fractions import Fraction
from typing import Dict, FrozenSet, Iterable, List, Optional, Sequence, Tuple

from .dvg import Digraph, QZ, Z


def rayleigh_difference(coeffs: Dict[FrozenSet[int], float], i: int, j: int, point: Dict[int, float]) -> float:
    """Delta_ij = d_i P d_j P - P d_i d_j P for P = sum_X coeffs[X] x^X, at a real point
    (coordinates i and j of the point are irrelevant)."""
    def part(a: bool, b: bool) -> float:
        s = 0.0
        for X, c in coeffs.items():
            if (i in X) != a or (j in X) != b:
                continue
            m = c
            for t in X:
                if t not in (i, j):
                    m *= point.get(t, 0.0)
            s += m
        return s
    return part(True, False) * part(False, True) - part(False, False) * part(True, True)


def lambda_of_path(D: Digraph, path: Sequence[int], exits: Optional[Dict[int, float]] = None):
    """lambda(P) = sum (d^+(t) + e_t) - back(P) (Lemma 9.7; e_t is the z^{-1} coefficient of the
    self-energy of t, 0 by default)."""
    e = exits or {}
    return sum(len(D.out[t]) + e.get(t, 0) for t in path) - D.back_arcs(path)


def lambda_exact(D: Digraph, path: Sequence[int], eta=None):
    """Read lambda off the exact Gamma(P) in Q(z): gamma = z^{|P|} Gamma(P) = 1 + lambda z^{-2} + ..."""
    import sympy as sp
    R = D.resolvent_function(eta=eta)
    g = D.gamma(path, R) * Z ** len(path)
    s = sp.Symbol('s')
    expr = sp.sympify(g.as_expr()).subs(sp.Symbol('z'), 1 / s)
    ser = sp.series(expr, s, 0, 3).removeO()
    c0, c1, c2 = (ser.coeff(s, k) for k in range(3))
    if c0 != 1 or c1 != 0:
        raise AssertionError('unexpected expansion %s' % ser)
    return c2


def shortest_path(D: Digraph, src: int, dst: int) -> Optional[List[int]]:
    prev = {src: None}
    queue = [src]
    while queue:
        u = queue.pop(0)
        if u == dst:
            break
        for w in D.out[u]:
            if w not in prev:
                prev[w] = u
                queue.append(w)
    if dst not in prev:
        return None
    p = [dst]
    while p[-1] != src:
        p.append(prev[p[-1]])
    return p[::-1]


def obstruction_certificate(D: Digraph, exact: bool = False):
    """For a digraph that is not SCC-symmetric, return (x, y, path, kappa_empty, kappa_R):
    x -> y is a one-way arc on a directed cycle, path = (y = q_0, ..., q_k = x) a shortest return
    path, and kappa computed from lambda of path sets (exactly from Gamma(P) if exact=True).
    Returns None for SCC-symmetric digraphs."""
    oneway = D.one_way_arcs_in_components()
    if not oneway:
        return None
    x, y = oneway[0]
    path = shortest_path(D, y, x)
    lam = (lambda p: lambda_exact(D, p)) if exact else (lambda p: lambda_of_path(D, p))
    k_empty = lam([x]) + lam([y]) - lam([x, y])
    inner = path[1:-1]
    k_R = lam(path[1:]) + lam(path[:-1]) - lam(path) - (lam(inner) if inner else 0)
    return x, y, path, k_empty, k_R


def first_order_factorization_holds(k: int) -> bool:
    """Symbolic check (sympy) that for G = prod_{t}(1 + v_t) + eps * E with E multi-affine with
    generic coefficients on k variables, d/d eps Delta_{01}(G) at eps = 0 equals Q * F."""
    import sympy as sp
    v = sp.symbols('v0:%d' % k)
    subsets = [frozenset(c) for r in range(k + 1) for c in itertools.combinations(range(k), r)]
    lam = {X: sp.Symbol('l_' + ('_'.join(map(str, sorted(X))) if X else 'e')) for X in subsets}
    eps = sp.Symbol('eps')
    G = sp.prod([1 + vi for vi in v]) + eps * sum(lam[X] * sp.prod([v[t] for t in X]) for X in subsets)
    D = sp.diff(G, v[0]) * sp.diff(G, v[1]) - G * sp.diff(G, v[0], v[1])
    first = sp.expand(sp.diff(D, eps).subs(eps, 0))
    R = range(2, k)
    Q = sp.prod([1 + v[r] for r in R])
    F = 0
    for r in range(len(R) + 1):
        for c in itertools.combinations(R, r):
            S = frozenset(c)
            F += (lam[S | {0}] + lam[S | {1}] - lam[S | {0, 1}] - lam[S]) * sp.prod([v[t] for t in S])
    return sp.expand(D.subs(eps, 0)) == 0 and sp.expand(first - Q * F) == 0


def triangle_band_value(zval) -> Fraction:
    """For the directed triangle, G(v) = sum_T gamma_T v^T and the pair (x, y) = (2, 0):
    Delta_{20}(v_1 = -1/2) = (4 - z^2) / (4 (z^2 - 2)^2), negative exactly for |z| > 2."""
    zv = Fraction(zval)
    s = zv * zv - 2
    g1 = (s + 1) / s
    g2 = (s + 2) / s
    v = Fraction(-1, 2)
    return (g1 * g1 - g2) + v * g2 * (g1 - 1) + v * v * g2 * (g2 - g1)
