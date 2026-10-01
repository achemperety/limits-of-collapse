#!/usr/bin/env bash
# Reproduces the computations of Sections 11.6-11.9. Same requirements as ../run_all.sh (passagemath, networkx,
# sympy). Pass --full to include the 20-minute sampled search of Proposition 11.15 and the 54-minute search for
# the 6-cycle with three frozen vertices; their logs from the original run are in logs/.
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

# 11.6 Two tails
run ge_verify       python3 ge_verify.py                  # Proposition 11.15: the four conflicts, independently
if [ "$1" = "--full" ]; then
  run ge_search     python3 ge_barrier.py 40000 7 11      # Proposition 11.15: pairs of 11.2 separated; sampled search
fi

# 11.7 Real environments for long cycles
run jensen          python3 meanfield_jensen.py           # Theorem 11.18: g3 and its non-real roots

# 11.8 Smallest environments
run resultants      python3 meanfield_complete.py         # Proposition 11.19: resultants for m = 3..10
run minenv_5_3      python3 min_env_indep.py 5 3          # Proposition 11.20
run minenv_6_1      python3 min_env_indep.py 6 1          # Proposition 11.20: the 6-cycle with U[C] empty
run minenv_6_2      python3 min_env_indep.py 6 2
if [ "$1" = "--full" ]; then
  run minenv_6_3    python3 min_env_indep.py 6 3
fi

# 11.9 Pruned cores of any size
run blowup_family   python3 blowup_check_family.py        # Theorem 11.22, sufficiency for m = 3..12
run twins           python3 twin_quadratic.py             # Theorem 11.22, Step 3: the twin equation is linear
run twin_steps      python3 twin_step3.py                 # Theorem 11.22, Steps 2-4: the products Gamma, symbolically
run blowup_zero     python3 blowup_zero.py                # Theorem 11.23 and Conjecture 11.27(1)
run cacti           python3 cactus_check.py               # Theorem 11.24 and Corollary 11.25
run census6         python3 census.py nonsym6_backarc.txt # Proposition 11.26(1)
run zero6           python3 zero7.py cores6.txt           # Proposition 11.26(2), six vertices
run zero7           python3 zero7.py cores7.txt           # Proposition 11.26(2), seven vertices
run cores7          python3 describe_cores.py cores7.txt  # the 107 seven-vertex pruned cores
