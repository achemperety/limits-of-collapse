#!/usr/bin/env bash
# Reproduces every row of the "Computational verification" table of the addendum (Section 10).
# Requirements: Python 3 with passagemath (module sage.all__sagemath_singular), networkx, numpy; sympy optional.
# Logs go to logs/; each step prints its wall time and the last lines of its log.
set -e
cd "$(dirname "$0")"
mkdir -p logs

run() {
  local name="$1"; shift
  local t0=$(date +%s)
  "$@" > "logs/$name.log" 2>&1
  echo "$name: $(( $(date +%s) - t0 )) s"
  tail -n 4 "logs/$name.log" | sed 's/^/    /'
}

run godsil        python3 godsil_check.py            # Lemma 10.C with indeterminate weights
run n3            python3 run_n3.py                  # 3 vertices, algebraic closure
run n4            python3 run_n4.py                  # 4 vertices, algebraic closure; writes n4_details.pkl
run corollary_3p  python3 check_C1.py                # Corollary 10.3' against the potentials found on 4 vertices
run gen5          python3 gen5.py                    # the 9,608 digraphs on 5 vertices; writes digraphs5.pkl
run cocycle       python3 valuation_cocycle.py       # valuation and set potentials on 3-4 vertices
run census        python3 setpot_symbolic.py         # set-potential census, 3-5 vertices, confirmed over Q(z)
run n5_complex    python3 search5_setpot.py          # Theorem 10.8, weights in the algebraic closure

t0=$(date +%s)                                       # Theorem 10.8, independent real-weight run on two cores
python3 search5.py 0 2 > logs/n5_real_0.log 2>&1 &
python3 search5.py 1 2 > logs/n5_real_1.log 2>&1
wait
echo "n5_real (two processes): $(( $(date +%s) - t0 )) s"
tail -n 1 logs/n5_real_0.log logs/n5_real_1.log | sed 's/^/    /'

run frozen        python3 env_test.py                # Proposition 10.10: frozen environments, 1 and 2 frozen vertices
run frozen_degs   python3 env_rational.py            # Proposition 10.10: degrees of the weights
run frozen_sextic python3 env_sextic.py              # Proposition 10.10: the sextic and its sum-of-squares form
run c4            python3 c4_proof_check2.py         # Proposition 10.7'': the three Rayleigh differences
run k1_classes    python3 oneway3.py                 # one arc: outcome classes over matching numbers
run k1_rule_a     python3 oneway_rule2.py 7 1200 12  # Theorem 10.11 on every position with a live arc
run k1_rule_b     python3 oneway_rule2.py 8 400 14
run k2_classes    python3 twoway.py 2 6000           # two arcs: conflicting classes
run k2_example    python3 k2_example.py              # Proposition 10.12
run ladders       python3 resdeg.py                  # Theorem 10.14
