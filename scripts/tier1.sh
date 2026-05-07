#!/usr/bin/env bash
# Tier-1 correctness check end-to-end:
#   1. Generate variant XMLs from examples/6_2_simulation_datastreams.xml
#   2. Build the project (clean) and compile the test/benchmark sources
#   3. Run mascotdatastreams.benchmark.Tier1Runner over all variants
#   4. Plot convergence + per-call wall time
#
# Outputs all land under sandbox/.
#
# Re-run this whole script any time after making code changes.

set -euo pipefail

cd "$(dirname "$0")/.."

PYTHON="${PYTHON:-conda run --no-capture-output -n biopython_env python}"

echo "==> [1/3] Generating variant XMLs"
$PYTHON scripts/tier1_generate_xmls.py

echo
echo "==> [2/3] Building + running Tier1Runner (ant clean tier1)"
ant clean tier1

echo
echo "==> [3/3] Plotting"
$PYTHON scripts/tier1_plot.py

echo
echo "Done. Outputs:"
echo "  sandbox/tier1/                 (variant XMLs)"
echo "  sandbox/tier1_results.csv      (raw logP per XML)"
echo "  sandbox/tier1_summary.txt      (text summary)"
echo "  sandbox/tier1_convergence.png  (plot)"
