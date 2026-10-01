"""Twin identity: for m = 2r and s > t >= 1, K(c(s,t)) = z K(d(s,t)), where c, d are the coupling
sequences of the fresh plays from an even and an odd class of C_m(s,t,s,t,...).  Verified as a
polynomial identity in (s, z) with s symbolic (valid for every integer s > t)."""
import sympy as sp, sys
s, z = sp.symbols('s z')
def couplings(t, r, start):
    m = 2 * r
    sz = lambda i: s if i % 2 == 0 else t
    b = []
    # game from class 'start' ends after t*m (+1 if start even) levels when s > t
    L = t * m + (1 if start == 0 else 0)
    for k in range(L - 1):
        b.append(sz(start + k + 1) - (k + 1) // m)
    return b
def K(b):
    # continuant of a path with len(b)+1 levels, all weights z, edge couplings b
    A, B = sp.Integer(1), z          # K of empty-after and single level
    # backward: K_k = z K_{k+1} - b_k K_{k+2}
    Kn2, Kn1 = sp.Integer(1), z
    for bk in reversed(b):
        Kn2, Kn1 = Kn1, sp.expand(z * Kn1 - bk * Kn2)
    return Kn1
ok = True
for t in range(1, 7):
    for r in range(2, 7):
        lhs = K(couplings(t, r, 0)); rhs = sp.expand(z * K(couplings(t, r, 1)))
        good = sp.expand(lhs - rhs) == 0
        ok &= good
        if not good: print('FAIL t=%d r=%d' % (t, r))
print('twin identity K(c) = z K(d) for all t<=6, m=2r<=12, symbolic s:', ok)
