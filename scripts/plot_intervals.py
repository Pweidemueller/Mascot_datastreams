#!/usr/bin/env python
"""
Plot the distribution of integration interval lengths from CSVs produced by
mascotdatastreams.benchmark.IntervalDumper.

For each dataset (small + SARS), reads the OLD-mode CSV and NEW-mode CSV,
sanity-checks that the tree-event positions (coalescent + sampling t_start
values) match between the two, then plots three panels:

  1. Tree-event-only interval lengths   (intervals between consecutive
     coalescent/sampling events; ignores grid + maxInterval subdivisions).
     Identical for OLD and NEW.
  2. OLD-path interval lengths          (tree events + mascotshifts grid).
  3. NEW-path interval lengths          (tree events + maxInterval splits).

Vertical lines mark the user-grid step and the maxInterval setting.

Inputs:  sandbox/intervals/<dataset>_{old,new}.csv
Outputs: sandbox/intervals_<dataset>.png

Usage:  conda run -n biopython_env python scripts/plot_intervals.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
INTERVALS = REPO / "sandbox" / "intervals"


def tree_event_times(df: pd.DataFrame) -> np.ndarray:
    """All distinct times at which a coalescent or sampling event happens.

    Includes 0.0 (the most-recent-tip side) plus every t_start where
    event_at_start is 'sampling' or 'coalescent'. Also picks up the very
    last tree event from the last row's t_end (since IntervalDumper records
    the event AT t_start, the final coalescence at the root only appears as
    the t_end of the last row, never as a t_start).
    """
    mask = df["event_at_start"].isin(["sampling", "coalescent"])
    times = df.loc[mask, "t_start"].to_numpy()
    times = np.concatenate([[0.0], times, [df["t_end"].iloc[-1]]])
    return np.unique(times)


def tree_only_interval_lengths(df: pd.DataFrame) -> np.ndarray:
    """Durations between consecutive tree events. Length = #tree_events - 1."""
    t = tree_event_times(df)
    diffs = np.diff(t)
    return diffs[diffs > 0]  # strip duplicate-time pairs


def all_interval_lengths(df: pd.DataFrame) -> np.ndarray:
    """All non-zero ODE step lengths (tree + grid for old; tree + maxInterval for new)."""
    return df.loc[df["length"] > 0, "length"].to_numpy()


def grid_step_estimate(df: pd.DataFrame) -> float:
    """Estimate the user's grid step from rows whose start is a grid breakpoint."""
    sub = df[df["event_at_start"] == "interval"]
    if len(sub) == 0:
        return float("nan")
    return float(sub["length"].median())


def plot_one(dataset_label: str, csv_old: Path, csv_new: Path,
             grid_step: float, max_interval: float, png_path: Path) -> None:
    df_old = pd.read_csv(csv_old)
    df_new = pd.read_csv(csv_new)

    # Sanity check: tree event times match between old and new for the same XML.
    t_old = tree_event_times(df_old)
    t_new = tree_event_times(df_new)
    same_count = len(t_old) == len(t_new)
    same_values = same_count and np.allclose(t_old, t_new, atol=1e-9)
    print(f"[{dataset_label}] tree events: old={len(t_old)} new={len(t_new)}  "
          f"{'OK' if same_values else 'MISMATCH'}")
    if not same_values:
        # Report the first mismatch for diagnosis but still proceed with the plot.
        n = min(len(t_old), len(t_new))
        for i in range(n):
            if abs(t_old[i] - t_new[i]) > 1e-9:
                print(f"  first mismatch at index {i}: old={t_old[i]} new={t_new[i]}")
                break

    tree_only = tree_only_interval_lengths(df_old)
    old_all = all_interval_lengths(df_old)
    new_all = all_interval_lengths(df_new)

    # Empirical grid step from the OLD CSV's "interval" rows (sanity-check vs.
    # the value the user passed in).
    empirical = grid_step_estimate(df_old)
    print(f"  empirical grid-step (median of OLD interval rows): {empirical:.6g} "
          f"(supplied: {grid_step:g})")
    print(f"  tree-only intervals: n={len(tree_only)}, "
          f"min={tree_only.min():.6g}, median={np.median(tree_only):.6g}, "
          f"max={tree_only.max():.6g}")
    print(f"  OLD all intervals:   n={len(old_all)}, "
          f"min={old_all.min():.6g}, median={np.median(old_all):.6g}, "
          f"max={old_all.max():.6g}")
    print(f"  NEW all intervals:   n={len(new_all)}, "
          f"min={new_all.min():.6g}, median={np.median(new_all):.6g}, "
          f"max={new_all.max():.6g}")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    # Shared log-x bin edges for fair visual comparison across panels.
    all_lengths = np.concatenate([tree_only, old_all, new_all])
    lo = max(all_lengths.min(), 1e-12)
    hi = all_lengths.max()
    bins = np.logspace(np.log10(lo), np.log10(hi), 60)

    def add_grid_lines(ax: plt.Axes) -> None:
        ax.axvline(grid_step, color="tab:red", linestyle="--", linewidth=1.4,
                   label=f"grid step ({grid_step:g})")
        # Slight offset only if maxInterval coincides with grid_step, so the
        # vertical lines don't completely hide each other.
        if abs(max_interval - grid_step) / max(grid_step, 1e-12) < 0.05:
            ax.axvline(max_interval, color="tab:blue", linestyle=":", linewidth=1.4,
                       label=f"maxInterval ({max_interval:g}) — coincides with grid")
        else:
            ax.axvline(max_interval, color="tab:blue", linestyle=":", linewidth=1.4,
                       label=f"maxInterval ({max_interval:g})")
        ax.set_xscale("log")
        ax.set_xlabel("interval length")
        ax.set_ylabel("count")
        ax.legend(loc="best", fontsize=8)
        ax.grid(True, which="both", alpha=0.3)

    ax = axes[0]
    ax.hist(tree_only, bins=bins, color="tab:gray", alpha=0.75)
    ax.set_title(f"{dataset_label}\ntree-event-only intervals (n={len(tree_only)})")
    add_grid_lines(ax)

    ax = axes[1]
    ax.hist(old_all, bins=bins, color="tab:red", alpha=0.65)
    ax.set_title(f"OLD path: tree + grid (n={len(old_all)})")
    add_grid_lines(ax)

    ax = axes[2]
    ax.hist(new_all, bins=bins, color="tab:blue", alpha=0.65)
    ax.set_title(f"NEW path: tree + maxInterval (n={len(new_all)})")
    add_grid_lines(ax)

    fig.tight_layout()
    fig.savefig(png_path, dpi=140)
    plt.close(fig)
    print(f"  wrote {png_path.relative_to(REPO)}")


def main() -> None:
    if not INTERVALS.is_dir():
        sys.exit(f"[plot_intervals] missing {INTERVALS}; run the dumper first")

    # SARS: mascotshifts step is 0.005, we matched maxInterval = 0.005.
    plot_one(
        dataset_label="SARS-CoV-2 1000seq",
        csv_old=INTERVALS / "sarscov2_old.csv",
        csv_new=INTERVALS / "sarscov2_new.csv",
        grid_step=0.005,
        max_interval=0.005,
        png_path=REPO / "sandbox" / "intervals_sarscov2.png",
    )

    # 6_2 small: gridRateShifts step is ~2e-4 (irregular; alternates 1e-4/2e-4).
    # We matched maxInterval = 2e-4.
    plot_one(
        dataset_label="6_2 small (125 seqs)",
        csv_old=INTERVALS / "6_2_small_old.csv",
        csv_new=INTERVALS / "6_2_small_new.csv",
        grid_step=2e-4,
        max_interval=2e-4,
        png_path=REPO / "sandbox" / "intervals_6_2_small.png",
    )


if __name__ == "__main__":
    main()
