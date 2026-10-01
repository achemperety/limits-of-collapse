# Verification scripts for §11, "Two arcs, the complex rescue, and the anatomy of set potentials"

**Run everything:** `./run_section11.sh` (about 16 minutes on two cores). Add `--full` to include the search of
Proposition 11.6 for the 5-cycle with two frozen vertices, which takes about 40 more minutes; its log from the
original run is `logs/minenv_5_2.log`. Requirements are those of `../run_all.sh` plus mpmath and sympy.

| Claim in §11 | Script | Expected output | Time |
| --- | --- | --- | --- |
| Theorem 11.1 (common tail), every position with the tail unvisited | `star_rule.py 1 800 11 3`, `2 1500 12 4`, `3 600 13 5` | 218,028 + 423,163 + 201,434 = 842,625 positions, 0 mismatches | 3 min |
| The six conflicting classes by geometry | `k2_conflicts.py` | 26,633 starts, 4,220 classes, conflicts: 4 chained, 2 common head | 4 min |
| Refinement by adjacency among key points | `k2_features.py 79 6000` | 12,107 classes, 3 conflicts | 1 min |
| Proposition 11.2 | `k2_barrier.py` | three pairs with identical key-local data and opposite outcomes | 1 s |
| Theorem 11.4 (the triangle), exact | `triangle_exact.py` | all identities; the resultant is the sextic; the sum-of-squares form | 1 s |
| Theorem 11.5 (every cycle) | `frozen_cycle.py` | F_m(alpha) = 0 for m = 3..8; residuals below 1e-55 on all branches for m = 3..7; exponents near -(m-2)/m; no real branch at z = 1e4, 1e8 | 1 min |
| Proposition 11.6 | `min_env.py 4 1`, `4 2`, `5 1` (and `5 2` with `--full`) | 276, 4,496, 6,560 (and 216,288) graphs U, none with a frozen potential | 1 min (+40) |
| Inside the band | `inband_scan.py`, `inband_intervals.py` | real branches: m = 3 whole band, m = 4 on (0.178, 2), m = 5 on (-1.024, 1.025), none for m = 6, 7 | 3 min |
| Lemma 11.8 against the exact test; census rows | `taxonomy.py` | 0 mismatches on all 8,172 non-SCC-symmetric digraphs on 3-5 vertices; 1 + 6 + 79 with a set potential | 30 s |
| Corollary 11.13 | `census_table.py` | the table of Corollary 11.13, total 79 | 1 s |
| Theorems 11.10 and 11.11, and three predictions | `families_check.py` | cycles m = 3..8 and books r = 2..6 as stated; the doubled 4-cycle realized by exits | 1 s |
| Formulas in the proof of Theorem 11.11 | `book_formulas.py` | True | 1 s |
| Theorem 11.12 | `eta_variety.py 3`, `eta_variety.py 4`, `core5_filter.py`, `eta_variety5.py` | the cores and self-energies of the table; 14 of 5,027 five-vertex cores pass the back-arc test | 3 min |
| Structure theorems on larger digraphs | `predictions.py` | 17 of 17 predictions confirmed | 1 s |
| Exit realizations of the exotic cores | `realizability.py` | 8,844 exit self-energies; bowtie 0, doubled triangle 0, 2-page book 22 | 5 s |

Shared modules: `dvg.py`, `setpot_symbolic.py`, `valuation_cocycle.py`, `oneway3.py` (random graphs), `env_test.py`
(frozen-environment test) and `gen5.py` (the 5-vertex digraphs).
