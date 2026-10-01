"""Blow-ups of directed cycles and the continuant criterion for set potentials (Section 13).

C_m(s_0, ..., s_{m-1}) replaces the i-th vertex of the directed m-cycle by a class S_i of s_i twins,
with an arc from every vertex of S_i to every vertex of S_{i+1} (indices mod m).

Theorem 13.9 (blow-up criterion).  With a constant self-energy e and y = z - e, a directed path that
starts in class j and has n vertices has Gamma(P) = K^(j)_n / K^(j)_0, where K^(j) are the backward
continuants of the branching sequence b^(j)_k = s_{j+k+1} - floor((k+1)/m).  Hence C_m(s) has a set
potential iff K^(j)_{cm} / K^(j)_0 does not depend on j for 1 <= c <= min(s).

Lemma 13.10 (twin identity) and Corollary 13.11: every two-periodic blow-up has a set potential.
"""
from __future__ import annotations

from typing import List, Sequence

from .dvg import Digraph, QZ, Z


def blowup_digraph(sizes: Sequence[int]) -> Digraph:
    """The digraph C_m(s_0, ..., s_{m-1}); class i occupies a contiguous block of vertices."""
    m = len(sizes)
    if m < 3 or min(sizes) < 1:
        raise ValueError('need m >= 3 classes of positive size')
    start = [0]
    for s in sizes:
        start.append(start[-1] + s)
    arcs = []
    for i in range(m):
        j = (i + 1) % m
        for a in range(start[i], start[i + 1]):
            for b in range(start[j], start[j + 1]):
                arcs.append((a, b))
    return Digraph(start[-1], arcs)


def branching(sizes: Sequence[int], j: int) -> List[int]:
    """b^(j)_k = s_{j+k+1} - floor((k+1)/m) for k = 0, 1, ..., up to and including the first zero."""
    m = len(sizes)
    b = []
    k = 0
    while True:
        bk = sizes[(j + k + 1) % m] - (k + 1) // m
        b.append(bk)
        if bk <= 0:
            b[-1] = 0
            return b
        k += 1


def backward_continuants(b: Sequence[int], e=None):
    """K_L = 1, K_{L+1} = 0, K_k = (z - e) K_{k+1} - b_k K_{k+2}; returns [K_0, ..., K_{L+1}] in Q(z)."""
    y = Z - (e if e is not None else QZ(0))
    L = len(b)
    K = [None] * (L + 2)
    K[L] = QZ(1)
    K[L + 1] = QZ(0)
    for k in range(L - 1, -1, -1):
        K[k] = y * K[k + 1] - b[k] * K[k + 2]
    return K


def gamma_formula(sizes: Sequence[int], j: int, n: int, e=None):
    """Gamma of any directed path with n vertices starting in class j (Theorem 13.9(1))."""
    K = backward_continuants(branching(sizes, j), e)
    return K[n] / K[0]


def criterion(sizes: Sequence[int], e=None) -> bool:
    """Theorem 13.9(2): set potential iff K^(j)_{cm}/K^(j)_0 is independent of j for c <= min(s)."""
    m = len(sizes)
    for c in range(1, min(sizes) + 1):
        vals = [gamma_formula(sizes, j, c * m, e) for j in range(m)]
        if any(v != vals[0] for v in vals[1:]):
            return False
    return True


def two_periodic(sizes: Sequence[int]) -> bool:
    m = len(sizes)
    return all(sizes[i] == sizes[(i + 2) % m] for i in range(m))


def twin_identity_holds(s: int, t: int, m: int, e=None) -> bool:
    """Lemma 13.10: for m even and s > t >= 1, K^(0)_k = (z - e) K^(1)_k for every even k <= t m."""
    if m % 2 or not s > t >= 1:
        raise ValueError('need m even and s > t >= 1')
    sizes = [s if i % 2 == 0 else t for i in range(m)]
    K0 = backward_continuants(branching(sizes, 0), e)
    K1 = backward_continuants(branching(sizes, 1), e)
    y = Z - (e if e is not None else QZ(0))
    return all(K0[k] == y * K1[k] for k in range(0, t * m + 1, 2))
