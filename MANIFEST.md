# Manifest

Every file in the archive, what it contains, and which results of the paper it supports. Theorem numbers refer to `paper/main_siam.pdf`.

## Paper (`paper/`)

| File | Contents |
| --- | --- |
| `main_siam.tex`, `preamble_siam.tex` | SIAM build (`siamart0516`; the class file is not included, see README) |
| `main_amsart.tex`, `preamble.tex` | `amsart` build, with identical numbering |
| `sections/front_siam.tex`, `sections/front_amsart.tex`, `sections/abstract.tex` | title, abstract, keywords, AMS classification |
| `sections/s00_intro.tex` | overview figure, contributions, status of results |
| `sections/s01.tex`–`s08.tex` | §§1–8: setting, local characterization, classical barriers, extensions, spectral formula, matching potentials (Theorem 6.10), algebraic collapse, exact boundary |
| `sections/s09.tex` | §9: real stability forces reciprocity (Theorem 9.5, Theorem 9.12, Proposition 9.14) |
| `sections/s10.tex` | §10: potentials over arbitrary fields (field independence, obstructions, Theorem 10.15) |
| `sections/s11.tex` | §11: the complex frontier; the counterexample (Theorem 11.6, Corollary 11.7, Propositions 11.8–11.10) |
| `sections/s12.tex` | §12: one-way arcs; star-tail rule (Theorem 12.2), barrier, $p$ tails |
| `sections/s13.tex` | §13: set potentials; cycles, books, pruned cores, blow-ups (Theorem 13.9, Lemma 13.10, Corollary 13.11), cacti (Theorem 13.13) |
| `sections/s14.tex` | §14: potentials versus bounded width (Theorems 14.1, 14.2, 14.6; Proposition 14.5) |
| `sections/s15.tex` | §15: status and open problems |
| `sections/app.tex` | Appendices: computational verification, notation, concordance with the earlier drafts |
| `sections/bib.tex` | bibliography |
| `limits_of_collapse_siam.tex`, `limits_of_collapse_amsart.tex` | the same paper as single self-contained files (all sections inlined) |
| `main_siam.pdf`, `main_amsart.pdf` | compiled builds (76 and 70 pages) |
| `Makefile` | `make amsart`, `make siam SIAMDIR=...` |

## Referee report (`review/`)

| File | Contents |
| --- | --- |
| `referee_report.tex`, `.pdf`, `.md` | Phase I referee report: 22 logical vulnerabilities, 7 missed opportunities, 7 formalization deficits, each mapped to its resolution; outcome of the research mandates |
| `make_refs.py`, `paper_refs.tex`, `paper_refs.json`, `paper.aux` | resolve the report's references from the paper's `.aux` file |
| `to_markdown.py` | Markdown export of the report |

## Python package (`src/dvgcollapse/`)

| Module | Contents | Paper |
| --- | --- | --- |
| `dvg.py` | `Digraph`: strong components (Tarjan), SCC-symmetry, reachable positions, outcomes, exact root resolvents over ℚ(z) (with self-energies), path products Γ, back arcs | §§5–6, Lemma 9.7 |
| `matching.py` | matching number, essential vertices, Gallai–Edmonds classes, weighted matching polynomials μ<sub>a</sub> by the vertex recurrence | §6 |
| `potentials.py` | the potential of Theorem 6.6; exact verification of a potential; set and valuation potentials | Thm 6.6, §13 |
| `real_stability.py` | Rayleigh differences, first-order data ℓ(P), exact and numerical obstruction certificates, the triangle band | §9 |
| `groebner_saturation.py` | the polynomial system of a (frozen) potential; saturated Gröbner test over ℚ(z) with Singular, or pure sympy with a Rabinowitsch variable; searches over graphs U | Lemma 10.1, Thm 10.15, Props 11.2, 11.10, 11.12–11.14 |
| `complex_counterexample.py` | D(m, r), Y<sub>m</sub>, prescribed elementary symmetric functions; `verify_symbolic` (every reachable position, exact) and `verify_numeric` (roots of Y<sub>m</sub>, 50 digits) | Thm 11.6, Props 11.8–11.9 |
| `star_tail.py` | the star-tail rule, exhaustive evaluator, fresh-entry outcomes, DVG solver for common tails, key-local data | Thm 12.2, Cors 12.3–12.4, Props 12.5–12.9 |
| `blowups.py` | blow-up digraphs, branching numbers, backward continuants, the continuant criterion, the twin identity | Thm 13.9, Lemma 13.10, Cor 13.11 |
| `treewidth.py` | ladders and closed ladders, resolvent degree, canonical subgames, the prime-tail comb, component-state potentials | Thms 14.1–14.6, Prop 14.5 |
| `solver.py` | unified exact DVG solver (common tails → component states → exhaustive); command line `dvg-solve` | Cors 6.8, 12.4, 14.7 |

## Tests (`tests/`, 41 tests)

| File | What it checks |
| --- | --- |
| `test_dvg.py` | positions and resolvents on cycles; strong components; Theorem 5.2 on random digraphs; Γ and back arcs |
| `test_potentials.py` | Theorem 6.6 on random SCC-symmetric digraphs; hierarchy examples (Theorem 13.1) |
| `test_groebner_saturation.py` | the three-vertex census; the rescued triangle; one frozen vertex never rescues it |
| `test_real_stability.py` | first-order factorization; certificates on all small non-SCC-symmetric digraphs; exact certificate; triangle band |
| `test_complex_counterexample.py` | D<sub>m</sub> is not SCC-symmetric; exact identities for (m, r) ∈ {(3,0), (4,0), (3,1)}; the form of Y<sub>3</sub>; residuals and non-real counts at large z |
| `test_star_tail.py` | the star rule against exhaustive search at every position of random instances; the DVG solver; the barrier examples of Propositions 12.5, 12.6 and 12.9 |
| `test_blowups.py` | Theorem 13.9(1) against all paths; the criterion against direct tests on all 115 blow-ups with ≤ 9 vertices; the twin identity; two-periodic families |
| `test_treewidth.py` | ladder degrees; component-state potentials; resolvents depend only on canonical subgames; the comb; closed ladders |
| `test_solver.py` | the unified solver against exhaustive search, fallback, command line |

## Scripts (`scripts/`; `run_all.sh` runs them)

| Script | Result | Needs |
| --- | --- | --- |
| `check_factorization.py`, `check_kappa.py` | Lemma 9.9; κ(∅) = 0, κ(R) = 1 | sympy |
| `check_c4.py` | the three Rayleigh expansions of Corollary 10.14 | sympy |
| `check_rescued_triangle.py` | Theorem 11.3 and the sextic of Proposition 11.2 | sympy |
| `complex_numeric_D3.py` | Theorem 11.6 for m = 3 at real and complex z | mpmath |
| `ym_arithmetic.py` | Proposition 11.9(3) | passagemath-pari |
| `search6_c3k3.py`, `search6_c3k3_rigorous.py` | Proposition 11.10 (1,928 graphs U, saturated over ℚ(z)) | passagemath-singular |
| `blowup_scan.py`, `blowup_identity.py`, `blowup_formula.py`, `twin_identity_symbolic.py` | Theorem 13.9, Lemma 13.10, Conjecture 13.16(1) on 331 size vectors | sympy |
| `comb_memory.py` | Proposition 14.5(2) for m ≤ 6 | sympy |
| `closed_ladders.py` | the closed-ladder table of §14.3 | sympy |
| `adversarial_*.py`, `ctrl_exact.py` | numerical sanity checks for §9 (not used in proofs): no real-stable completion of cycle path data at large z; a symmetric control is found | numpy, scipy |
| `run_all.out` | output of the last run | |

## Legacy verification (`verification/legacy/`)

The verification package of the earlier drafts (`dvg_addendum_verification-3`), unchanged: `run_all.sh` and the folders `section11/` and `section11b/` reproduce every row of the addendum's verification tables (Theorems 10.8, 10.11, 11.1–11.27 in the old numbering; see Appendix C of the paper for the new numbers). It needs passagemath-singular.
