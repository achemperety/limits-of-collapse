# Verification scripts for §10, "Matching potentials and one-way arcs: an addendum"

These scripts reproduce every row of the addendum's "Computational verification" table. The scripts for §11.1–11.4
(two arcs, the complex rescue, set potentials) are in the folder `section11/`, and those for §11.6–11.9 (two tails,
real environments, smallest environments, pruned cores of any size) in `section11b/`, each with its own README and
runner.

**Requirements.** Python 3 with passagemath's Singular interface (the module `sage.all__sagemath_singular`),
networkx and numpy. sympy is optional; it is used only for a cross-check in `env_sextic.py`.

**Run everything.** `./run_all.sh` writes one log per step to `logs/` and prints each step's wall time. The
complete run takes about 17 minutes on a two-core machine (measured times below).

**Exactness.** All algebra is exact: rational arithmetic, or Gröbner bases over the field ℚ(z) with saturation,
which decide solvability over the algebraic closure of ℚ(z). Evaluation at z = 37/11 is used only to refute an
identity. That is rigorous because every resolvent denominator is a monic integer polynomial and 37/11 is not an
algebraic integer, so it is never a pole.

## Library

| File | Contents |
| --- | --- |
| `dvg.py` | digraphs, SCC-symmetry, reachable positions, exact resolvents R(W, v) over ℚ(z), weighted matching polynomials μₐ(U[S]) |
| `solver.py` | `has_potential(n, arcs, U)`: the ideal of μₐ(U[W − v]) − R(W, v)·μₐ(U[W]), saturated by every μₐ(U[W]); unit ideal ⇔ no potential |
| `valuation_cocycle.py` | valuation potentials, and set potentials tested at z = 37/11 (refutation only) |
| `setpot_symbolic.py` | `set_potential_exact(n, arcs)`: set potentials decided exactly over ℚ(z) |

## What each script checks

| Claim in §10 | Command | Expected output | Time |
| --- | --- | --- | --- |
| Lemma 10.C with indeterminate weights | `python3 godsil_check.py` | `tests 1680 failures 0` | 1 s |
| Potentials on 3 vertices (algebraic closure) | `python3 run_n3.py` | 13 SCC-symmetric digraphs have one; the 3 others have none | 2 s |
| Potentials on 4 vertices (algebraic closure) | `python3 run_n4.py` | `{(True, True): 104, (False, False): 114}`; writes `n4_details.pkl` | 45 s |
| Corollary 10.3′ against the 310 potentials on 4 vertices | `python3 check_C1.py` | 389 incomparable pairs, 198 symmetric pairs inside a component, 0 violations | < 1 s |
| The 9,608 digraphs on 5 vertices | `python3 gen5.py` | `9608`; writes `digraphs5.pkl` | 6 s |
| Set-potential census (Theorem 10.13) | `python3 setpot_symbolic.py` | 86 flagged (1, 6, 79 on 3, 4, 5 vertices), all 86 confirmed over ℚ(z) | 6 s |
| Theorem 10.8, weights in the algebraic closure | `python3 search5_setpot.py` | 8,055 non-SCC-symmetric; 79 with a set potential; 80,798 pairs (D, U); 0 potentials | 6 min |
| Theorem 10.8, independent real-weight run | `python3 search5.py 0 2` and `python3 search5.py 1 2` in parallel | 183,442 + 183,109 = 366,551 pairs; 19,805 + 19,574 = 39,379 pass the exact filter; 0 potentials | 2 min on two cores |
| Proposition 10.10, frozen environments | `python3 env_test.py` | directed triangle: 0 graphs U with one frozen vertex, 4 with two; the other two triples: 0 and 0 | 10 s |
| Proposition 10.10, degrees of the weights | `python3 env_rational.py` | degree 3 on the cycle, 6 on the environment, for all four U | 1 s |
| Proposition 10.10, the sextic | `python3 env_sextic.py` | both environment weights have the sextic as minimal polynomial; sum-of-squares identity holds | 1 s |
| Proposition 10.7″, the three Rayleigh differences | `python3 c4_proof_check2.py` | Δ₁₂, Δ₁₃, Δ₀₂ as stated in the proof | 1 s |
| One arc: outcome classes over matching numbers | `python3 oneway3.py` | `instances 58412 classes 42 conflicts 0` | 3.5 min |
| Theorem 10.11 on every position with a live arc | `python3 oneway_rule2.py 7 1200 12` and `python3 oneway_rule2.py 8 400 14` | 204,142 + 68,257 = 272,399 positions, 0 mismatches | 1 min |
| Two arcs: conflicting classes (Proposition 10.12) | `python3 twoway.py 2 6000` | `k 2 instances 26633 classes 4220 conflicts 6` | 3.5 min |
| Proposition 10.12, the explicit example | `python3 k2_example.py` | identical deficit vectors, different outcomes | < 1 s |
| Theorem 10.14, ladders | `python3 resdeg.py` | same-direction rails: degrees 6, 14, 30, 62, 126, 254 for m = 2, …, 7 | 1 s |
