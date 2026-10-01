import mpmath as mp
import frozen_cycle as fc
mp.mp.dps = 30
def allreal(m, zr, k):
    th = (mp.acos(zr / 2) + 2 * mp.pi * k) / m
    tc = 2 * mp.cos(th)
    al = sum(mp.mpf(int(c)) * tc ** i for i, c in enumerate(fc.cheb_U(m).list()))
    coeffs = [sum(mp.binomial(r, j) * (-al) ** (r - j) * fc.muP(j, zr) for j in range(r + 1)) / mp.factorial(r) for r in range(m)]
    # all roots real <=> Sturm-free check via numerical roots
    b = mp.polyroots(list(reversed(coeffs)), maxsteps=300, extraprec=150)
    return all(abs(mp.im(x)) < mp.mpf(10) ** -15 for x in b)
for m in (4, 5):
    N = 4000; cur = None; segs = []
    for i in range(1, N):
        zr = -2 + 4 * mp.mpf(i) / N
        ks = tuple(k for k in range(m) if allreal(m, zr, k))
        if ks != cur:
            segs.append([float(zr), ks]); cur = ks
    print('m = %d: (left endpoint, branches with all-real weights)' % m)
    for s in segs: print('   z >= %.4f : %s' % (s[0], s[1]))
