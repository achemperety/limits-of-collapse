"""Decide existence of a local weighted matching potential (U, a) for a digraph D, rigorously over
K = Q(z): the unknown weights alpha_x range over the algebraic closure of K. The system is
   mu_a(U[W-v]) - R(W,v) mu_a(U[W]) = 0   for every reachable position (W, v),
with mu_a(U[W]) != 0 imposed by saturation. Unit ideal <=> no potential (for this U), even with
algebraic-function weights; in particular none with weights in Q(z)."""
import itertools, sys, time
from sage.all__sagemath_singular import PolynomialRing, QQ
import dvg

def potential_ideal(n, arcs, Uedges, field=None, zval=None):
    K = dvg.Kz if field is None else field
    A = PolynomialRing(K, ['a%d' % i for i in range(n)], order='degrevlex')
    al = A.gens()
    P = dvg.reachable_positions(n, arcs)
    R = dvg.make_resolvent(n, arcs, field=K, zval=zval)
    memo = {}
    eqs = []; Ws = set()
    for (W, v) in P:
        num = dvg.mu_weighted(W - {v}, Uedges, al, A, memo)
        den = dvg.mu_weighted(W, Uedges, al, A, memo)
        eqs.append(num - A(R(W, v)) * den)
        Ws.add(W)
    I = A.ideal([e for e in eqs if e != 0])
    nonzero = [dvg.mu_weighted(W, Uedges, al, A, memo) for W in Ws]
    return A, I, nonzero

def has_potential(n, arcs, Uedges, field=None, zval=None, return_ideal=False):
    A, I, nonzero = potential_ideal(n, arcs, Uedges, field, zval)
    if I.is_zero():
        return (True, I) if return_ideal else True
    J = I
    for f in nonzero:
        if f.is_constant():
            if f == 0: return (False, None) if return_ideal else False
            continue
        J = J.saturation(A.ideal([f]))[0]
        if J.is_one():
            return (False, J) if return_ideal else False
    ok = not J.is_one()
    return (ok, J) if return_ideal else ok

def all_graphs_on(n):
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
    for k in range(len(pairs) + 1):
        for c in itertools.combinations(pairs, k):
            yield frozenset(c)
