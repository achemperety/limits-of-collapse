#!/usr/bin/env bash
# Reproduce the computations of the revision.  Quick checks run by default (a few minutes);
# pass --full for the long ones (Groebner search over 1,928 graphs, adversarial searches).
set -euo pipefail
cd "$(dirname "$0")"
PY=${PYTHON:-python3}
run() { echo "=== $*"; "$@"; }
run $PY check_factorization.py
run $PY check_kappa.py
run $PY check_c4.py
run $PY check_rescued_triangle.py
run $PY complex_numeric_D3.py
run $PY blowup_scan.py
run $PY twin_identity_symbolic.py
run $PY comb_memory.py 5
run $PY closed_ladders.py 6
if $PY -c "import sage.all__sagemath_pari" 2>/dev/null; then run $PY ym_arithmetic.py; else echo "skip ym_arithmetic.py (needs passagemath-pari)"; fi
if [[ "${1:-}" == "--full" ]]; then
  if $PY -c "import sage.all__sagemath_singular" 2>/dev/null; then
    run $PY search6_c3k3.py
    run $PY search6_c3k3_rigorous.py
  else
    echo "skip search6_c3k3*.py (needs passagemath-singular)"
  fi
  run $PY adversarial_cycle.py 5 8 0 800
  run $PY adversarial_de.py 5 8 1 800 0
  run $PY ctrl_exact.py
fi
echo "done"
