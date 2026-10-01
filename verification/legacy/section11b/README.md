# Verification scripts for §11.6–11.9 of the addendum

**Run everything:** `./run_section11b.sh` (about 5 minutes on one core). Add `--full` to include the sampled
search of Proposition 11.15 and the search for the 6-cycle with three frozen vertices, about 75 more minutes; their
logs from the original run are in `logs/`. Requirements are those of `../run_all.sh` (passagemath with Singular,
networkx, sympy). The census inputs were generated with nauty (`geng`, `directg`) and the two C filters here; the
generation commands are listed at the end and are not rerun by the script.

| Claim | Script | Expected output |
| --- | --- | --- |
| Proposition 11.15(2): the four conflicts, checked independently (brute-force matchings and outcomes) | `ge_verify.py` | four pairs with the same GE-local data and outcomes True, False |
| Proposition 11.15(1) and the sampled search | `ge_barrier.py 40000 7 11` (`--full`) | the pairs of Proposition 11.2 have different GE classes; 38,805 two-tail starts, 13,680 classes, 4 conflicts |
| Theorem 11.18: the Jensen polynomial g3 and its roots | `meanfield_jensen.py` | closed form True; root identity True; two non-real roots at z = 3 for random real alpha |
| Proposition 11.19: resultants for symmetric environments | `meanfield_complete.py` | nonzero for m = 3..10 with U[C] empty (±w^a (w^2 − 1)^b) and complete |
| Proposition 11.20: the 5-cycle, three frozen vertices, U[C] empty (4 minutes) | `min_env_indep.py 5 3` | 9,168 orbits, none with a frozen potential |
| Proposition 11.20: the 6-cycle with U[C] empty | `min_env_indep.py 6 1`, `6 2`, `6 3` (`--full`) | 14, 724 and 59,976 orbits, none with a frozen potential (the last takes 54 minutes) |
| Theorem 11.22: sufficiency for m = 3..12 | `blowup_check_family.py` | all rotation continuants equal for every m |
| Theorem 11.22, Step 3: the twin equation is linear | `twin_quadratic.py` | degree 1 for (m, r) = (3, 2), (3, 3), (4, 2), (5, 2) |
| Theorem 11.22, Steps 2-4: the products Gamma and the twin equation | `twin_step3.py` | True in all 6 + 7 cases (m up to 7, r up to 3), all weights symbolic |
| Theorem 11.23 and Conjecture 11.27(1): blow-ups with zero self-energies | `blowup_zero.py` | True exactly for the two-periodic size sequences: 12 of the 20 blow-ups tried (the list repeats four of them in rotated form) |
| Theorem 11.24 and Corollary 11.25: cacti | `cactus_check.py` | predicted: True; perturbed: False; all eta = 0: False, for all 8 cacti |
| Proposition 11.26(1): the six-vertex census | `census.py nonsym6_backarc.txt` | 2,291 of 9,113 with a set potential; the table of pruned cores |
| Proposition 11.26(2): zero self-energies on six and seven vertices | `zero7.py cores6.txt`, `zero7.py cores7.txt` | 5 of 79 and 3 of 469; C6, C3(2,2,2), C4(1,2,1,2) and C7 among the pruned cores |
| The 107 seven-vertex pruned cores | `describe_cores.py cores7.txt` | 15 blow-ups of cycles and 92 other cores |

Shared modules: `dvg.py` (games, resolvents over Q(z)), `setpot_symbolic.py` (the exact set-potential test),
`setpot_eta.py` (the same test with self-energies), `valuation_cocycle.py`, `env_test.py` (frozen environments),
`describe_cores.py` (pruned cores and blow-up recognition), `census.py` (components and pruning), `oneway3.py`
(random graphs for the sampled search).

Generation of the census inputs (nauty 2.8, about 107 minutes for seven vertices):

```
gcc -O2 -o backarc backarc.c && gcc -O2 -o backarc2 backarc2.c
geng -q 6 | directg -T -q | ./backarc2 > nonsym6_backarc.txt   # not SCC-symmetric, back-arc test passed
geng -c -q 6 | directg -T -q | ./backarc > cores6.txt           # strongly connected, back-arc test passed
geng -c -q 7 | directg -T -q | ./backarc > cores7.txt
```
