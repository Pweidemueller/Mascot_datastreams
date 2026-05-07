#!/usr/bin/env bash
# Tier-2 timing benchmark end-to-end:
#   1. Generate variant XMLs at two grid baselines (user1x and user5x) plus the
#      multi-parameter perturbation set.
#   2. Build, compile-test, run Tier2Runner over the main sweep.
#   3. Run Tier2Runner over the multi-parameter perturbations.
#   4. Plot per-(dataset, baseline) panels and the multi-param bias plot.
#
# Outputs:
#   sandbox/tier2/                       (main-sweep variant XMLs)
#   sandbox/tier2_results.csv            (main-sweep raw timings)
#   sandbox/tier2_<dataset>_<base>.png   (per-baseline panels)
#   sandbox/tier2_summary.txt            (text summary)
#   sandbox/tier2_multiparam/            (perturbation XMLs)
#   sandbox/tier2_multiparam_results.csv
#   sandbox/tier2_multiparam.png
#
# Tunables:
#   TIMED=N    timed iterations per XML in main sweep (default 30)
#   WARMUP=N   warm-up iterations per XML in main sweep (default 5)

set -euo pipefail

cd "$(dirname "$0")/.."

PYTHON="${PYTHON:-conda run --no-capture-output -n biopython_env python}"
WARMUP="${WARMUP:-5}"
TIMED="${TIMED:-30}"

echo "==> [1/5] Generating main-sweep XMLs"
$PYTHON scripts/tier2_generate_xmls.py

echo
echo "==> [2/5] Generating multi-parameter perturbation XMLs"
$PYTHON scripts/tier2_multiparam.py

echo
echo "==> [3/5] Building + running main sweep (warmup=$WARMUP, timed=$TIMED)"
ant clean tier2 -Dtier2.warmup=$WARMUP -Dtier2.timed=$TIMED

echo
echo "==> [4/5] Running multi-parameter sweep (no rebuild)"
# tier2-multiparam-only skips the dep chain so install-dependencies isn't
# re-invoked (which would otherwise fail because deps/ is already populated).
ant tier2-multiparam-only

echo
echo "==> [5/5] Plotting"
$PYTHON scripts/tier2_plot.py
$PYTHON scripts/tier2_multiparam.py --plot

echo
echo "Done."
