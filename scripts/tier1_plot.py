#!/usr/bin/env python
"""
Plot Tier-1 convergence: logP from the new (max-interval) path at decreasing
maxInterval, against the old path at the original and densified grids.

Reads:  sandbox/tier1_results.csv
Writes: sandbox/tier1_convergence.png
        sandbox/tier1_summary.txt

CSV schema (written by Tier1Runner.java):
  xml_path,mode,parameter,logP,wall_ms

Usage:  python scripts/tier1_plot.py
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parent.parent
CSV_PATH = REPO / "sandbox" / "tier1_results.csv"
PNG_PATH = REPO / "sandbox" / "tier1_convergence.png"
TXT_PATH = REPO / "sandbox" / "tier1_summary.txt"


def parse_max_interval(param: str) -> float | None:
    if not param.startswith("max"):
        return None
    return float(param[3:])


def main() -> None:
    if not CSV_PATH.is_file():
        sys.exit(f"[tier1_plot] missing {CSV_PATH} — run `ant tier1` first")

    rows = list(csv.DictReader(CSV_PATH.open()))
    if not rows:
        sys.exit("[tier1_plot] CSV is empty")

    # Split rows by mode
    old_grid1x = next((r for r in rows if r["mode"] == "old" and r["parameter"] == "grid1x"), None)
    old_grid5x = next((r for r in rows if r["mode"] == "old" and r["parameter"] == "grid5x"), None)
    new_rows = sorted(
        [r for r in rows if r["mode"] == "new"],
        key=lambda r: parse_max_interval(r["parameter"]) or 0.0,
        reverse=True,
    )

    if not new_rows:
        sys.exit("[tier1_plot] no new-mode rows in CSV; nothing to plot")

    xs = [parse_max_interval(r["parameter"]) for r in new_rows]
    ys = [float(r["logP"]) for r in new_rows]
    walls = [float(r["wall_ms"]) for r in new_rows]

    ref_old = float(old_grid5x["logP"]) if old_grid5x else None
    ref_orig = float(old_grid1x["logP"]) if old_grid1x else None

    # Convergence plot: logP vs maxInterval, log x-axis
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    ax = axes[0]
    ax.plot(xs, ys, marker="o", label="new (useMaxInterval=true)")
    if ref_old is not None:
        ax.axhline(ref_old, color="tab:red", linestyle="--", label=f"old grid (5x dense): {ref_old:.4f}")
    if ref_orig is not None:
        ax.axhline(ref_orig, color="tab:gray", linestyle=":", label=f"old grid (original): {ref_orig:.4f}")
    ax.set_xscale("log")
    ax.invert_xaxis()
    ax.set_xlabel("maxInterval (smaller → more accurate)")
    ax.set_ylabel("logP at fixed parameters")
    ax.set_title("Tier 1: logP convergence")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, which="both", alpha=0.3)

    # Wall-time plot
    ax = axes[1]
    ax.plot(xs, walls, marker="s", color="tab:orange", label="new path")
    if old_grid1x is not None:
        ax.axhline(float(old_grid1x["wall_ms"]), color="tab:gray", linestyle=":",
                   label=f"old grid (1x): {float(old_grid1x['wall_ms']):.2f} ms")
    if old_grid5x is not None:
        ax.axhline(float(old_grid5x["wall_ms"]), color="tab:red", linestyle="--",
                   label=f"old grid (5x): {float(old_grid5x['wall_ms']):.2f} ms")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.invert_xaxis()
    ax.set_xlabel("maxInterval")
    ax.set_ylabel("wall time per logP call (ms)")
    ax.set_title("Tier 1: per-call wall time")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    PNG_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(PNG_PATH, dpi=140)
    print(f"[tier1_plot] wrote {PNG_PATH.relative_to(REPO)}")

    # Text summary
    lines = ["Tier 1 — logP convergence summary", "=" * 38, ""]
    if ref_old is not None:
        lines.append(f"Reference (old grid, 5x dense):   {ref_old:.6f}")
    if ref_orig is not None:
        lines.append(f"Old grid (original):              {ref_orig:.6f}")
        if ref_old is not None:
            lines.append(f"  |old(1x) - old(5x)|           = {abs(ref_orig - ref_old):.6f}")
    lines.append("")
    lines.append("New path (useMaxInterval=true):")
    for r in new_rows:
        mi = parse_max_interval(r["parameter"])
        lp = float(r["logP"])
        diff = (lp - ref_old) if ref_old is not None else float("nan")
        lines.append(f"  maxInterval={mi:>9.5g}  logP={lp:.6f}  diff_vs_old5x={diff:+.6f}  wall_ms={float(r['wall_ms']):.2f}")
    TXT_PATH.write_text("\n".join(lines) + "\n")
    print(f"[tier1_plot] wrote {TXT_PATH.relative_to(REPO)}")
    print()
    print("\n".join(lines))


if __name__ == "__main__":
    main()
