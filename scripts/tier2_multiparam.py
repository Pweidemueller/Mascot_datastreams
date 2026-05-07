#!/usr/bin/env python
"""
Multi-parameter bias-stability spot check.

For each dataset, generate small perturbations of `migrationRatesSkyline`
around its initial value (1.0). For each perturbed XML, produce two variants:

  * old   — useMaxInterval=false (user's published grid)
  * new   — useMaxInterval=true at the maxInterval that approximately matches
            the user's mascotshifts step (drop-in target).

Run all of them via Tier2Runner. Plot logP_new − logP_old across the
perturbations: a roughly constant offset means the new path's bias tracks
the old path's bias as parameters change, i.e. they should produce
indistinguishable MCMC posteriors. A varying offset means the methods bias
the posterior differently.

The SARS XML's `mascotshifts` step is 5e-3 → drop-in maxInterval = 5e-3.
The 6_2 small XML's `gridRateShifts` step is 2e-4 → drop-in maxInterval = 2e-4.

Usage:  conda run -n biopython_env python scripts/tier2_multiparam.py
        ant clean tier2-multiparam
        conda run -n biopython_env python scripts/tier2_multiparam.py --plot
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUT_DIR = REPO / "sandbox" / "tier2_multiparam"
CSV_PATH = REPO / "sandbox" / "tier2_multiparam_results.csv"
PNG_PATH = REPO / "sandbox" / "tier2_multiparam.png"

DATASETS: list[tuple[str, Path, float]] = [
    # (label, src_xml, drop-in maxInterval matching the user's mascotshifts step)
    ("6_2_small",
     REPO / "examples" / "6_2_simulation_datastreams.xml",
     2e-4),
    ("sarscov2_1000seq",
     REPO / "examples" /
     "SARSCoV2_Epsilon_BayArea_results_1000seq_datastreams_noclip.xml",
     5e-3),
]

PERTURBATIONS = [
    ("mig0.5", 0.5),
    ("mig1.0", 1.0),
    ("mig2.0", 2.0),
    ("mig5.0", 5.0),
    ("mig0.1", 0.1),
]

MASCOT_PATTERN = re.compile(
    r'(<distribution\s+id="Mascot\.t:SimDataset"\s+spec="mascotdatastreams\.distribution\.MascotLogPflag"[^>]*?)(>)',
)
COUPLED_RUN_PATTERN = re.compile(
    r'<run\s+id="mcmc"\s+spec="coupledMCMC\.CoupledMCMC"([^>]*?)>',
)
MIG_PATTERN = re.compile(
    r'(<parameter\s+id="migrationRatesSkyline\.t:SimDataset"[^>]*>)([^<]*)(</parameter>)',
)


def patch_run_spec(t: str) -> str:
    def repl(m: re.Match[str]) -> str:
        cl = re.search(r'chainLength="([^"]+)"', m.group(1))
        return f'<run id="mcmc" spec="MCMC" chainLength="{cl.group(1) if cl else "1"}">'
    return COUPLED_RUN_PATTERN.sub(repl, t)


def patch_mascot(t: str, use_max_interval: bool, maxInterval: float | None) -> str:
    m = MASCOT_PATTERN.search(t)
    if not m:
        sys.exit("[tier2_multiparam] no MascotLogPflag tag")
    extra = ""
    if use_max_interval:
        extra = f' implementation="java" useMaxInterval="true" maxInterval="{maxInterval}"'
    return t[: m.start()] + m.group(1) + extra + m.group(2) + t[m.end():]


def patch_migration(t: str, scale: float) -> str:
    """Multiply the body of the migrationRatesSkyline parameter by `scale`."""
    m = MIG_PATTERN.search(t)
    if not m:
        sys.exit("[tier2_multiparam] no migrationRatesSkyline parameter")
    open_tag, body, close_tag = m.group(1), m.group(2), m.group(3)
    values = body.split()
    new_values = [f"{float(v) * scale:g}" for v in values] if len(values) > 1 \
        else [f"{float(values[0]) * scale:g}"]
    return t[: m.start()] + open_tag + " ".join(new_values) + close_tag + t[m.end():]


def generate() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    n = 0
    for dataset, src, mi_match in DATASETS:
        if not src.is_file():
            print(f"[tier2_multiparam] missing {src}, skipping {dataset}")
            continue
        base = patch_run_spec(src.read_text())
        for pert_label, scale in PERTURBATIONS:
            scaled = patch_migration(base, scale)
            for mode_label, use_mi, mi in [
                ("old", False, None),
                ("new", True, mi_match),
            ]:
                text = patch_mascot(scaled, use_mi, mi)
                out = OUT_DIR / f"{dataset}__{pert_label}__{mode_label}__match.xml"
                out.write_text(text)
                n += 1
                print(f"  wrote {out.relative_to(REPO)}")
    print(f"[tier2_multiparam] {n} variants written")


def plot() -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if not CSV_PATH.is_file():
        sys.exit(f"[tier2_multiparam] missing {CSV_PATH} — run the harness first")
    rows = list(csv.DictReader(CSV_PATH.open()))

    by_pair: dict[tuple[str, str], dict[str, dict]] = defaultdict(dict)
    for r in rows:
        # CSV columns: dataset, baseline=<pert_label>, mode (old/new), parameter (=match)
        key = (r["dataset"], r["baseline"])
        by_pair[key][r["mode"]] = r

    # Group rows by dataset to plot per-dataset panels.
    per_dataset: dict[str, list[tuple[str, float, float, float]]] = defaultdict(list)
    for (dataset, pert_label), pair in sorted(by_pair.items()):
        if "old" not in pair or "new" not in pair:
            continue
        old_lp = float(pair["old"]["logP"])
        new_lp = float(pair["new"]["logP"])
        diff = new_lp - old_lp
        per_dataset[dataset].append((pert_label, old_lp, new_lp, diff))

    fig, axes = plt.subplots(1, len(per_dataset), figsize=(6 * len(per_dataset), 4.6),
                             squeeze=False)
    for ax, (dataset, entries) in zip(axes[0], sorted(per_dataset.items())):
        # Ensure consistent x-order: as listed in PERTURBATIONS
        order = {p[0]: i for i, p in enumerate(PERTURBATIONS)}
        entries.sort(key=lambda e: order.get(e[0], 99))
        labels = [e[0] for e in entries]
        diffs = [e[3] for e in entries]
        xs = list(range(len(labels)))
        ax.plot(xs, diffs, marker="o", color="tab:purple")
        ax.set_xticks(xs)
        ax.set_xticklabels(labels, rotation=30, ha="right")
        ax.set_ylabel("logP(new) − logP(old)")
        ax.set_title(f"{dataset}\nbias stability across migration-rate scales")
        ax.axhline(0, color="black", linewidth=0.6, alpha=0.3)
        ax.grid(True, alpha=0.3)
        # Annotate logP values
        for x, e in zip(xs, entries):
            ax.annotate(f"old={e[1]:.2f}\nnew={e[2]:.2f}",
                        (x, e[3]), textcoords="offset points",
                        xytext=(0, 10), fontsize=7, ha="center")
    fig.tight_layout()
    fig.savefig(PNG_PATH, dpi=140)
    print(f"[tier2_multiparam] wrote {PNG_PATH.relative_to(REPO)}")

    # Text summary
    print()
    print("Multi-parameter bias check")
    print("=" * 30)
    for dataset, entries in sorted(per_dataset.items()):
        order = {p[0]: i for i, p in enumerate(PERTURBATIONS)}
        entries.sort(key=lambda e: order.get(e[0], 99))
        diffs = [e[3] for e in entries]
        spread = max(diffs) - min(diffs)
        mean = sum(diffs) / len(diffs)
        print(f"\n{dataset}:")
        print(f"  {'pert':<10} {'logP_old':>14} {'logP_new':>14} {'Δ=new-old':>14}")
        for label, old_lp, new_lp, diff in entries:
            print(f"  {label:<10} {old_lp:>14.4f} {new_lp:>14.4f} {diff:>+14.4f}")
        print(f"  mean Δ = {mean:.4f}, spread (max-min) = {spread:.4f}")
        print(f"  → bias is {'stable' if spread < 1.0 else 'NOT stable'}: spread {spread:.2f} logP across perturbations")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--plot", action="store_true",
                   help="If set, parse the CSV and plot. Otherwise generate XMLs.")
    args = p.parse_args()
    if args.plot:
        plot()
    else:
        generate()


if __name__ == "__main__":
    main()
