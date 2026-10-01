"""For which real z does some branch of the K_{m,m-1} family give all-real weights?"""
import mpmath as mp
from frozen_cycle import cheb_U, muP
mp.mp.dps = 40
def weights(m, zr, k):
    th = (mp.acos(zr / 2) + 2 * mp.pi * k) / m if abs(zr) < 2 else None
    tc = 2 * mp.cos(th)
    al = sum(mp.mpf(int(c)) * tc ** i for i, c in enumerate(cheb_U(m).list()))
    coeffs = [sum(mp.binomial(r, j) * (-al) ** (r - j) * muP(j, zr) for j in range(r + 1)) / mp.factorial(r) for r in range(m)]
    betas = mp.polyroots(list(reversed(coeffs)), maxsteps=300, extraprec=200)
    return tc, al, betas
for m in range(3, 8):
    N = 800
    good = []
    for i in range(1, N):
        zr = -2 + 4 * mp.mpf(i) / N
        for k in range(m):
            tc, al, betas = weights(m, zr, k)
            if all(abs(mp.im(b)) < 1e-25 for b in betas):
                good.append((float(zr), k))
    ks = sorted(set(k for _, k in good))
    frac = len(set(z for z, _ in good)) / (N - 1)
    print('m = %d: fraction of in-band z with an all-real branch: %.3f; branches used: %s' % (m, frac, ks))
