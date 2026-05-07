#!/usr/bin/env python
"""
Plot Tier-2 timing benchmark.

For each (dataset, baseline) we report:
  * logP (the model output)
  * n_doEuler_per_call (deterministic mechanism metric — primary)
  * wall_ms (secondary; varies with machine load)

The reference logP for the diff column is **the OLD path at the same baseline
(grid 1x in the user-published config)**. That's the practical drop-in target:
"can the new path produce the user's published-grid logP, faster?". A separate
"densest-config" diff is also reported for context.

Inputs:  sandbox/tier2_results.csv
Outputs: sandbox/tier2_<dataset>_<baseline>.png
         sandbox/tier2_summary.txt
"""
from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parent.parent
CSV_PATH = REPO / "sandbox" / "tier2_results.csv"
TXT_PATH = REPO / "sandbox" / "tier2_summary.txt"


def parse_max_interval(param: str) -> float | None:
    return float(param[3:]) if param.startswith("max") else None


def annotate_row(r: dict, logP_old: float, logP_dense: float) -> dict:
    return {
        "param": r["parameter"],
        "logP": float(r["logP"]),
        "diff_old": float(r["logP"]) - logP_old,
        "diff_dense": float(r["logP"]) - logP_dense,
        "n_euler": float(r["n_doEuler_per_call"]),
        "mean_ms": float(r["mean_ms"]),
        "std_ms": float(r["std_ms"]),
    }


def main() -> None:
    if not CSV_PATH.is_file():
        sys.exit(f"[tier2_plot] missing {CSV_PATH} — run `ant tier2` first")

    rows = list(csv.DictReader(CSV_PATH.open()))
    if not rows:
        sys.exit("[tier2_plot] CSV is empty")

    by_group: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for r in rows:
        by_group[(r["dataset"], r["baseline"])].append(r)

    summary_lines: list[str] = ["Tier 2 — timing summary", "=" * 38, ""]

    for (dataset, baseline), grp_rows in sorted(by_group.items()):
        old_row = next((r for r in grp_rows if r["mode"] == "old"), None)
        new_rows = sorted(
            [r for r in grp_rows if r["mode"] == "new"],
            key=lambda r: parse_max_interval(r["parameter"]) or 0.0,
            reverse=True,
        )
        if old_row is None:
            print(f"[tier2_plot] no old config in {dataset}/{baseline}, skipping")
            continue

        logP_old = float(old_row["logP"])
        # densest = config with most n_euler in this group
        densest = max(grp_rows, key=lambda r: float(r["n_doEuler_per_call"]))
        logP_dense = float(densest["logP"])

        old_a = annotate_row(old_row, logP_old, logP_dense)
        new_a = [(parse_max_interval(r["parameter"]),
                  annotate_row(r, logP_old, logP_dense)) for r in new_rows]

        # ---- Plot ----
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
        title = f"{dataset} | {baseline}"

        # Panel 1: logP convergence vs maxInterval. Solid = new path.
        # Reference (old grid) and densest config drawn as horizontal lines.
        ax = axes[0]
        if new_a:
            xs = [mi for mi, _ in new_a]
            ys = [a["logP"] for _, a in new_a]
            ax.plot(xs, ys, marker="o", color="tab:blue", label="new path")
        ax.axhline(logP_old, color="tab:red", linestyle="--",
                   label=f"old @ {baseline} ({logP_old:.3f})")
        if densest["mode"] != "old":
            ax.axhline(logP_dense, color="black", linestyle=":",
                       label=f"densest ({densest['mode']} {densest['parameter']}, {logP_dense:.3f})")
        ax.set_xscale("log")
        ax.invert_xaxis()
        ax.set_xlabel("maxInterval")
        ax.set_ylabel("logP")
        ax.set_title(f"{title} — logP vs maxInterval")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(loc="best", fontsize=7)

        # Panel 2: n_doEuler_per_call vs maxInterval.
        ax = axes[1]
        if new_a:
            xs = [mi for mi, _ in new_a]
            ys = [a["n_euler"] for _, a in new_a]
            ax.plot(xs, ys, marker="o", color="tab:blue", label="new path")
        ax.axhline(old_a["n_euler"], color="tab:red", linestyle="--",
                   label=f"old @ {baseline}: {old_a['n_euler']:.0f}")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.invert_xaxis()
        ax.set_xlabel("maxInterval")
        ax.set_ylabel("doEuler calls per likelihood eval")
        ax.set_title(f"{title} — ODE step count")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(loc="best", fontsize=7)

        # Panel 3: cost vs work — ms vs n_euler. Both methods overlaid.
        ax = axes[2]
        if new_a:
            xs = [a["n_euler"] for _, a in new_a]
            ys = [a["mean_ms"] for _, a in new_a]
            errs = [a["std_ms"] for _, a in new_a]
            ax.errorbar(xs, ys, yerr=errs, marker="o", linestyle="-",
                        color="tab:blue", label="new path")
            for x, y, (mi, _) in zip(xs, ys, new_a):
                ax.annotate(f"{mi:g}", (x, y), textcoords="offset points",
                            xytext=(4, 4), fontsize=7)
        ax.errorbar([old_a["n_euler"]], [old_a["mean_ms"]], yerr=[old_a["std_ms"]],
                    marker="s", linestyle="none", color="tab:red",
                    label=f"old @ {baseline} ({old_a['mean_ms']:.2f} ms)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("doEuler calls per likelihood eval")
        ax.set_ylabel("wall time per call (ms)")
        ax.set_title(f"{title} — cost vs work")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(loc="best", fontsize=8)

        fig.tight_layout()
        png = REPO / "sandbox" / f"tier2_{dataset}_{baseline}.png"
        fig.savefig(png, dpi=140)
        plt.close(fig)
        print(f"[tier2_plot] wrote {png.relative_to(REPO)}")

        # ---- Text summary block ----
        summary_lines.append(f"### {dataset}  baseline={baseline}")
        summary_lines.append("")
        summary_lines.append(
            f"  match target (old path at this baseline): logP = {logP_old:.4f}"
        )
        summary_lines.append(
            f"  densest config (best logP estimate):      logP = {logP_dense:.4f}"
            f"  ({densest['mode']} {densest['parameter']})"
        )
        summary_lines.append("")
        summary_lines.append(
            f"  {'config':<20} {'logP':>14} {'Δold':>9} {'Δdense':>9}"
            f" {'n_euler':>10} {'mean ms':>10} {'std ms':>9}"
        )
        for label, a in [(f"old @ {baseline}", old_a)] + [
            (f"new max={mi:g}", a) for mi, a in new_a
        ]:
            summary_lines.append(
                f"  {label:<20} {a['logP']:>14.4f} {a['diff_old']:>+9.3f}"
                f" {a['diff_dense']:>+9.3f} {a['n_euler']:>10.1f}"
                f" {a['mean_ms']:>10.3f} {a['std_ms']:>9.3f}"
            )

        # Matched-accuracy block: smallest |diff_old| among new configs.
        if new_a:
            best_match = min(new_a, key=lambda x: abs(x[1]["diff_old"]))
            mi, a = best_match
            ratio_n = a["n_euler"] / max(old_a["n_euler"], 1.0)
            ratio_ms = a["mean_ms"] / max(old_a["mean_ms"], 1.0)
            summary_lines.append("")
            summary_lines.append(
                f"  → drop-in match (smallest |Δold|): new max={mi:g} → Δold={a['diff_old']:+.3f}"
            )
            summary_lines.append(
                f"     n_euler ratio new/old = {ratio_n:.2f}×    wall ratio new/old = {ratio_ms:.2f}×"
                f"  ({'speedup' if ratio_ms < 1 else 'slowdown'} of {1/ratio_ms if ratio_ms<1 else ratio_ms:.2f}×)"
            )
        summary_lines.append("")

    TXT_PATH.write_text("\n".join(summary_lines) + "\n")
    print(f"[tier2_plot] wrote {TXT_PATH.relative_to(REPO)}")
    print()
    print("\n".join(summary_lines))


if __name__ == "__main__":
    main()
