"""Proposition 11.9(3): irreducibility, discriminant factors, and Galois groups of specializations of Y_m
(m = 3..6).  Requires the PARI interface of passagemath (sage.all__sagemath_pari)."""
import sympy as sp, time
from sage.all__sagemath_pari import pari
z, y = sp.symbols('z y')
def muP(k):
    a, b = sp.Integer(1), z
    if k == 0: return a
    for _ in range(k-1): a, b = b, sp.expand(z*b - a)
    return b
def Ym(m):
    # mu(P_m) * Y_m(y), cleared of the common factor z where possible
    expr = sp.expand(muP(m)*y**(m+1) - sp.factorial(m)*z*sum(muP(l)*y**l/sp.factorial(l) for l in range(m+1)))
    expr = sp.expand(sp.cancel(expr/sp.gcd_list(sp.Poly(expr, y).all_coeffs())))
    return expr
for m in range(3, 7):
    t0 = time.time()
    P = Ym(m)
    Pp = sp.Poly(P, y, z)
    irr = Pp.is_irreducible
    lc = sp.factor(sp.Poly(P, y).LC())
    disc = sp.Poly(sp.discriminant(P, y), z)
    fac = sp.factor_list(disc.as_expr())
    simple = [ (f, e) for f, e in fac[1] if e % 2 == 1 and sp.Poly(f, z).degree() > 0]
    # a simple irreducible factor f (exponent 1) that is squarefree and coprime to the leading coefficient
    ok_simple = [f for f, e in fac[1] if e == 1 and sp.Poly(f, z).degree() > 0 and sp.gcd(f, sp.Poly(P, y).LC()) == 1]
    gal = []
    for z0 in (3, 5, 7):
        spec = sp.Poly(sp.expand(P.subs(z, z0)), y)
        s = str(spec.as_expr()).replace('**', '^')
        gal.append(str(pari('polgalois(%s)' % s)))
    print('m=%d irreducible=%s  LC=%s  disc factors (deg,exp)=%s  simple factors=%d  Galois(z0=3,5,7)=%s  (%.1fs)' % (
        m, irr, lc, [(sp.Poly(f, z).degree(), e) for f, e in fac[1]], len(ok_simple), gal, time.time()-t0), flush=True)
