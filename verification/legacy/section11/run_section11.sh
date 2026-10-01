#!/usr/bin/env bash
# Reproduces every computation of Section 11. Same requirements as ../run_all.sh (passagemath, networkx, numpy,
# mpmath, sympy). Pass --full to include the 40-minute search of Proposition 11.6 for the 5-cycle with |F| = 2.
set -e
cd "$(dirname "$0")"
mkdir -p logs

run() {
  local name="$1"; shift
  local t0=$(date +%s)
  "$@" > "logs/$name.log" 2>&1
  echo "$name: $(( $(date +%s) - t0 )) s"
  tail -n 3 "logs/$name.log" | sed 's/^/    /'
}

run gen5            python3 gen5.py                   # the 9,608 digraphs on 5 vertices (digraphs5.pkl)

# 11.1 Two one-way arcs
run star_a          python3 star_rule.py 1 800 11 3   # Theorem 11.1 on every position
run star_b          python3 star_rule.py 2 1500 12 4
run star_c          python3 star_rule.py 3 600 13 5
run k2_conflicts    python3 k2_conflicts.py           # the six conflicting classes by geometry
run k2_adjacency    python3 k2_features.py 79 6000    # refinement by adjacency among key points
run k2_barrier      python3 k2_barrier.py             # Proposition 11.2

# 11.2 The complex rescue
run triangle        python3 triangle_exact.py         # Theorem 11.4, exact
run cycles          python3 frozen_cycle.py           # Theorem 11.5
run minenv_4_1      python3 min_env.py 4 1            # Proposition 11.6
run minenv_4_2      python3 min_env.py 4 2
run minenv_5_1      python3 min_env.py 5 1
if [ "$1" = "--full" ]; then run minenv_5_2 python3 min_env.py 5 2; fi
run band            python3 inband_scan.py            # inside the band: numerical scan
run band_intervals  python3 inband_intervals.py

# 11.3 Set potentials
run components      python3 taxonomy.py               # Lemma 11.8 against the exact test; census rows
run census          python3 census_table.py           # Corollary 11.13
run families        python3 families_check.py         # Theorems 11.10 and 11.11, and three predictions
run book_formulas   python3 book_formulas.py          # the explicit formulas in the proof of Theorem 11.11
run cores_3         python3 eta_variety.py 3          # Theorem 11.12
run cores_4         python3 eta_variety.py 4
run cores_5_filter  python3 core5_filter.py
run cores_5         python3 eta_variety5.py
run predictions     python3 predictions.py            # structure theorems on larger digraphs
run realizability   python3 realizability.py          # exit realizations of the exotic cores
