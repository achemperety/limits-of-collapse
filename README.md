# limits-of-collapse

Companion archive for the revised paper **“Alternation, Invariants, and the Limits of Collapse”** (Ruslan Baiandin, revised September 30, 2026). It contains the paper, in SIAM (`siamart0516`) and `amsart` builds, the Phase I referee report with a response map, a tested Python package, `dvgcollapse`, and every script behind the paper's computational claims.

## Main results of the revision and how they are established

| Result | Paper | Status | Where it is checked |
| --- | --- | --- | --- |
| A digraph has a local weighted matching potential over ℚ(z), ℝ(z) or real Puiseux series **iff** it is SCC-symmetric (formerly Conjecture 6.10) | Thm 6.10, Thm 9.5 | proof | `tests/test_real_stability.py`, `scripts/check_factorization.py`, `scripts/check_kappa.py` |
| Real frozen environments never rescue a strongly connected non-symmetric digraph (formerly Conjecture 10.9′) | Thm 9.12 | proof | same |
| The triangle's real band \|z\| < 2 is sharp for **every** environment | Prop 9.14 | proof | `tests/test_real_stability.py` |
| **Over the algebraic closure of ℚ(z) the statement is false**: the directed m-cycle plus m+1 isolated vertices, U = K<sub>m,m+1</sub>, zero cycle weights, and environment weights at the roots of Y<sub>m</sub> | Thm 11.6, Cor 11.7 | proof + exact and numerical checks | `dvgcollapse.complex_counterexample`, `tests/test_complex_counterexample.py`, `scripts/complex_numeric_D3.py` |
| Y<sub>m</sub> is unramified at ∞; exactly m−1 (odd m) or m (even m) weights are non-real at large real z; its monodromy group is S<sub>m+1</sub> for 3 ≤ m ≤ 6 | Prop 11.9 | proof; item (3) computer-assisted | `scripts/ym_arithmetic.py` |
| Triangle plus k isolated vertices has a complex potential iff k ≥ 4 | Prop 11.10 | computer-assisted (1,928 saturated ideals over ℚ(z)) | `scripts/search6_c3k3*.py` |
| DVG is polynomial when the one-way arcs of every strong component share a tail (the star-tail rule) | Thm 12.2, Cor 12.4 | proof | `dvgcollapse.star_tail`, `tests/test_star_tail.py` |
| Continuant criterion for blow-ups of directed cycles; twin identity; every two-periodic blow-up has a set potential | Thm 13.9, Lemma 13.10, Cor 13.11 | proof | `dvgcollapse.blowups`, `tests/test_blowups.py`, `scripts/blowup_scan.py` |
| Directed cacti (the cacti taxonomy) | Thm 13.13, Cor 13.14 | proof | `verification/legacy/section11b` |
| Exact resolvents: degree ≥ 2<sup>m</sup>−1 for every multiplicative potential however indexed; 2<sup>m</sup> distinct resolvents at one vertex at treewidth 3; linear-size component-state circuits | Thm 14.2, Prop 14.5, Thm 14.6 | proof | `dvgcollapse.treewidth`, `tests/test_treewidth.py`, `scripts/comb_memory.py`, `scripts/closed_ladders.py` |

Theorem numbers refer to the SIAM build (`paper/main_siam.pdf`); the amsart build has the same numbering.

## Layout

```
paper/            LaTeX sources (shared sections, two preambles), compiled PDFs, Makefile
review/           Phase I referee report (.tex, .pdf, .md) and the scripts that resolve its references
src/dvgcollapse/  the Python package
tests/            pytest suite (41 tests, about 45 s)
scripts/          one-off computations of the revision; run_all.sh reproduces them
verification/legacy/  the verification package of the earlier drafts, unchanged
MANIFEST.md       file-by-file manifest
```

## Install and test

```bash
python -m pip install -e ".[test,numeric]"
pytest -q                      # 41 tests
scripts/run_all.sh             # quick reproduction (a few minutes)
scripts/run_all.sh --full      # adds the Groebner search and the numerical adversarial runs
```

Exact Gröbner bases with saturation and PARI Galois groups need the optional passagemath modules (`pip install passagemath-singular passagemath-pari`). Scripts that need them skip themselves when the modules are missing. The `groebner_saturation` module falls back to a pure-sympy Rabinowitsch test when Singular is not installed.

## Using the package

```python
from dvgcollapse.dvg import Digraph
from dvgcollapse.solver import outcome
from dvgcollapse import complex_counterexample as cc

D = Digraph(5, [(0, 1), (1, 0), (1, 2), (2, 3), (3, 1), (3, 4)])
print(outcome(D, 0))                 # (True/False, method used)

print(cc.Y_cleared(3))               # the environment polynomial of D_3
print(cc.verify_symbolic(3))         # every potential identity of D_3, exactly
```

From the command line: `dvg-solve --n 3 --arcs "0-1 1-2 2-0" --start 0`.

## Building the paper

```bash
make -C paper amsart               # needs a standard TeX Live
make -C paper siam SIAMDIR=/path/to/siam/macros
```

The SIAM class `siamart0516.cls` is **not** included. Its license forbids distributing it apart from the complete SIAM macro package, which is available from SIAM and on CTAN (`macros/latex/contrib/siam`). The compiled PDFs of both builds are in `paper/`.
