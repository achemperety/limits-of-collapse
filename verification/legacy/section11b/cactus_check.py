"""Checks of the cactus theorem (directed cacti: every block a chordless directed cycle) with the exact test.

Positive cases use self-energies predicted by the theorem; each is followed by a perturbed copy, which must fail.
"""
from setpot_eta import setpot_eta, cycles_digraph, z, Kz


def mu_path(k, A):
    a, b = Kz(1), Kz(A)            # mu(P_0), mu(P_1)
    if k == 0:
        return a
    for _ in range(k - 1):
        a, b = b, A * b - a
    return b


def f(m, A):
    """root resolvent of the directed path on m - 1 vertices, all with weight A (a petal of length m entered after h)"""
    return mu_path(m - 2, A) / mu_path(m - 1, A)


cases = []
# 1. bowtie: two triangles sharing h = 0 (Theorem 11.12 table: eta(h) = e - 1/((z-e) - 1/(z-e)), e = 0)
n, arcs = cycles_digraph([[0, 1, 2], [0, 3, 4]])
cases.append(('bowtie', n, arcs, {0: -f(3, z)}))
# 2. flower of three triangles: eta(h) = -2 f_3(z)
n, arcs = cycles_digraph([[0, 1, 2], [0, 3, 4], [0, 5, 6]])
cases.append(('flower of 3 triangles', n, arcs, {0: -2 * f(3, z)}))
# 3. chain of three triangles, cut vertices 0 and 3: each cut vertex lies on two triangles
n, arcs = cycles_digraph([[0, 1, 2], [0, 3, 4], [3, 5, 6]])
cases.append(('chain of 3 triangles', n, arcs, {0: -f(3, z), 3: -f(3, z)}))
# 4. flower of two 4-cycles with genuinely period-2 petals: beta = z, gamma = z - 1/z
beta, gamma = z, z - 1 / z
RB = (gamma * beta - 1) / (beta * (gamma * beta - 2))
wh = gamma + RB
n, arcs = cycles_digraph([[0, 1, 2, 3], [0, 4, 5, 6]])
cases.append(('flower (4,4), period 2', n, arcs, {0: z - wh, 1: z - beta, 2: z - gamma, 3: z - beta,
                                                 4: z - beta, 5: z - gamma, 6: z - beta}))
# 5. flower of a triangle and a 4-cycle, explicit rational solution (alpha = z on the triangle)
alpha = z
cA = alpha - f(3, alpha)
beta = (z ** 2 + 1 - 1 / z ** 2) / cA
gamma = (z ** 2 + 2) / beta
wh = gamma + f(3, alpha)
RB = (gamma * beta - 1) / (beta * (gamma * beta - 2))
assert wh == alpha + RB
n, arcs = cycles_digraph([[0, 1, 2], [0, 3, 4, 5]])
cases.append(('flower (3,4)', n, arcs, {0: z - wh, 1: z - alpha, 2: z - alpha, 3: z - beta, 4: z - gamma, 5: z - beta}))
# 6. chain of three 4-cycles, uniform: cut vertices on two cycles get -f_4(z)
n, arcs = cycles_digraph([[0, 1, 2, 3], [0, 4, 5, 6], [5, 7, 8, 9]])
cases.append(('chain of three 4-cycles', n, arcs, {0: -f(4, z), 5: -f(4, z)}))
# 7. a 5-cycle with a 5-cycle attached at two of its vertices (uniform, 13 vertices)
n, arcs = cycles_digraph([[0, 1, 2, 3, 4], [0, 5, 6, 7, 8], [2, 9, 10, 11, 12]])
cases.append(('5-cycle with two 5-cycles attached', n, arcs, {0: -f(5, z), 2: -f(5, z)}))
# 8. flower of a triangle, a 4-cycle and a 5-cycle is not uniform; skip. A star of four triangles:
n, arcs = cycles_digraph([[0, 1, 2], [0, 3, 4], [0, 5, 6], [0, 7, 8]])
cases.append(('flower of 4 triangles', n, arcs, {0: -3 * f(3, z)}))

for name, n, arcs, eta in cases:
    ok = setpot_eta(n, arcs, eta)
    bad = dict(eta); v0 = max(bad, key=lambda v: v) if bad else 0
    bad[v0] = bad.get(v0, 0) + 1 / z ** 3
    ko = setpot_eta(n, arcs, bad)
    zero = setpot_eta(n, arcs, {})
    print('%-36s n=%2d  predicted: %s | perturbed: %s | all eta = 0: %s' % (name, n, ok, ko, zero), flush=True)
