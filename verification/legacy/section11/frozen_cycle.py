"""Frontier II. Frozen potentials for the directed m-cycle with a complete bipartite environment K_{m, m-1}.

Claim: with t a root of 2 T_m(t/2) = z and alpha = U_m(t/2) (the matching polynomial of the path P_m at t),
the environment weights are the m-1 roots of
    Q(x) = sum_{r < m} F_r(alpha) x^r / r!,   F_r(a) = sum_k C(r,k) (-a)^(r-k) mu(P_k)(z),
and (U = K_{C,F}, cycle weights alpha, environment weights = roots of Q) is a frozen potential.

Checks: (1) F_m(alpha) = 0 identically in t; (2) the frozen identity at every reachable position of C_m, evaluated
with the actual weighted matching polynomials of U at random complex z (mpmath, 60 digits); (3) Puiseux scaling
of the environment weights at large real z, and whether any branch is real."""
import itertools, random, cmath
import mpmath as mp
from sage.all__sagemath_singular import PolynomialRing, QQ, binomial

mp.mp.dps = 60
Pt = PolynomialRing(QQ, 't'); t = Pt.gen()


def cheb_T2(m):     # 2 T_m(t/2)
    a, b = Pt(2), t
    if m == 0: return a
    for _ in range(m - 1): a, b = b, t * b - a
    return b


def cheb_U(m):      # U_m(t/2) = mu(P_m)(t)
    a, b = Pt(1), t
    if m == 0: return a
    for _ in range(m - 1): a, b = b, t * b - a
    return b


def muP(k, zz):
    a, b = 1, zz
    if k == 0: return a
    for _ in range(k - 1): a, b = b, zz * b - a
    return b


def F(r, a, zz):
    return sum(binomial(r, k) * (-a) ** (r - k) * muP(k, zz) for k in range(r + 1))




def cheb_roots(m, zc):
    c = [mp.mpc(int(x)) for x in cheb_T2(m).list()]   # low -> high
    c[0] -= zc
    return mp.polyroots(list(reversed(c)), maxsteps=400, extraprec=400)

def mu_weighted(S, E, w):
    S = tuple(sorted(S))
    if not S: return mp.mpc(1)
    v, rest = S[0], S[1:]
    val = w[v] * mu_weighted(rest, E, w)
    for u in rest:
        if (v, u) in E: val -= mu_weighted(tuple(x for x in rest if x != u), E, w)
    return val


def check_potential(m, zc, tc):
    """all branches: tc a root of 2T_m(t/2) = zc."""
    al = sum(mp.mpc(c) * tc ** i for i, c in enumerate(cheb_U(m).list()))
    coeffs = []
    for r in range(m):
        s = mp.mpc(0)
        for k in range(r + 1):
            s += mp.binomial(r, k) * (-al) ** (r - k) * muP(k, zc)
        coeffs.append(s / mp.factorial(r))
    betas = mp.polyroots(list(reversed(coeffs)), maxsteps=200, extraprec=200)
    C = list(range(m)); Fv = list(range(m, 2 * m - 1))
    E = {(i, j) for i in C for j in Fv}
    w = {i: al for i in C}
    for j, b in zip(Fv, betas): w[j] = b
    # resolvents of the directed m-cycle game: interval of length L starting at v: mu(P_{L-1})/mu(P_L)
    worst = 0
    for start in C:
        for L in range(m, 0, -1):          # unvisited interval of length L starting at start
            W = [(start + i) % m for i in range(L)]
            R = muP(L - 1, zc) / muP(L, zc)
            lhs = mu_weighted(tuple(W[1:]) + tuple(Fv), E, w) / mu_weighted(tuple(W) + tuple(Fv), E, w)
            worst = max(worst, abs(lhs - R))
    return worst, al, betas




if __name__ == '__main__':
    print('(1) F_m(alpha) = 0 identically, with z = 2T_m(t/2), alpha = U_m(t/2):')
    for m in range(3, 9):
        zz, al = cheb_T2(m), cheb_U(m)
        print('   m = %d:' % m, F(m, al, zz) == 0)


    print('(2) frozen identity at every reachable position, random complex z, every branch t:')
    random.seed(3)
    for m in range(3, 8):
        zc = mp.mpc(random.uniform(2.5, 4), random.uniform(-1, 1))
        ts = cheb_roots(m, zc)
        errs = [check_potential(m, zc, tc)[0] for tc in ts]
        print('   m = %d: %d branches, max residual %s' % (m, len(ts), mp.nstr(max(errs), 3)))

    print('(3) large real z: scaling of the environment weights and reality')
    for m in range(3, 9):
        out = []
        for zr in (mp.mpf(10) ** 4, mp.mpf(10) ** 8):
            ts = cheb_roots(m, zr)
            branch_info = []
            for tc in ts:
                _, al, betas = check_potential(m, mp.mpc(zr), tc) if m <= 5 else (0, None, None)
                if betas is None:
                    al = sum(mp.mpc(c) * tc ** i for i, c in enumerate(cheb_U(m).list()))
                    coeffs = []
                    for r in range(m):
                        s = sum(mp.binomial(r, k) * (-al) ** (r - k) * muP(k, zr) for k in range(r + 1))
                        coeffs.append(s / mp.factorial(r))
                    betas = mp.polyroots(list(reversed(coeffs)), maxsteps=300, extraprec=300)
                allreal = abs(mp.im(al)) < 1e-20 and all(abs(mp.im(b)) < 1e-20 for b in betas)
                branch_info.append((allreal, max(abs(b) for b in betas)))
            out.append((zr, branch_info))
        (z1, b1), (z2, b2) = out
        expo = [mp.log(bb2 / bb1) / mp.log(z2 / z1) for (_, bb1), (_, bb2) in zip(b1, b2)]
        print('   m = %d: any branch with all weights real: %s | fitted exponent of max|beta|: %s (predicted %s)'
              % (m, any(r for r, _ in b1 + b2), mp.nstr(max(expo), 5), mp.nstr(-mp.mpf(m - 2) / m, 5)))
