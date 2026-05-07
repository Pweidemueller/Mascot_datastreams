#!/usr/bin/env python
"""
Compare the OLD and NEW Tier-3 MCMC runs.

Reads the BEAST .log files from sandbox/tier3_old/ and sandbox/tier3_new/,
plus the run-time hours/Msample numbers from the run.log files, and reports:
  * Wall-time speedup (samples per unit time).
  * Per-parameter posterior summary (mean, 95% HPD-ish quantiles) on the
    post-burn-in tail.
  * ESS estimates with a simple autocorrelation-based estimator.
  * Whether the two posteriors overlap meaningfully on the headline parameters.
"""
from __future__ import annotations

from pathlib import Path
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent

OLD_DIR = REPO / "sandbox" / "tier3_old"
NEW_DIR = REPO / "sandbox" / "tier3_new"
OLD_LOG = OLD_DIR / "SARSCoV2_Epsilon_BayArea_results_1000seq_datastreams_noclip.log"
NEW_LOG = NEW_DIR / "SARSCoV2_Epsilon_BayArea_results_1000seq_datastreams_noclip_interval.log"
OLD_RUN = OLD_DIR / "run.log"
NEW_RUN = NEW_DIR / "run.log"


def load_trace(p: Path) -> pd.DataFrame:
    return pd.read_csv(p, sep="\t", comment="#", engine="python")


def last_progress(run_log: Path) -> tuple[int, str]:
    """Return (last sample number, hours/Msamples string) from BEAST run.log."""
    last = None
    pat = re.compile(r"^\s*(\d+)\s.*\s(\d+h\d+m\d+s)/Msamples\s*$")
    with run_log.open() as f:
        for line in f:
            m = pat.match(line.rstrip())
            if m:
                last = (int(m.group(1)), m.group(2))
    return last if last else (0, "?")


def hms_to_hours(s: str) -> float:
    m = re.match(r"(\d+)h(\d+)m(\d+)s", s)
    if not m:
        return float("nan")
    h, mn, sec = (int(x) for x in m.groups())
    return h + mn / 60 + sec / 3600


def autocorr_ess(x: np.ndarray) -> float:
    """Effective sample size via initial-positive-sequence estimator (Geyer)."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    if n < 4:
        return float(n)
    x = x - x.mean()
    var = (x * x).mean()
    if var == 0:
        return float(n)
    # Compute autocorrelation up to where it goes negative.
    lags = []
    cum = 1.0
    for k in range(1, min(n - 1, 10000)):
        rho = (x[:-k] * x[k:]).mean() / var
        if rho <= 0 and k > 1:
            break
        lags.append(rho)
        cum += 2 * rho
    return float(n / max(cum, 1.0))


def summary(df: pd.DataFrame, params: list[str], burn_frac: float = 0.2) -> pd.DataFrame:
    """Per-parameter mean, 95% (2.5/97.5) quantiles, ESS, on post-burn-in samples."""
    n = len(df)
    burn = int(n * burn_frac)
    out = []
    for p in params:
        if p not in df.columns:
            continue
        x = df[p].iloc[burn:].to_numpy()
        out.append({
            "param": p,
            "mean": np.mean(x),
            "q025": np.quantile(x, 0.025),
            "q975": np.quantile(x, 0.975),
            "ess": autocorr_ess(x),
            "n": len(x),
        })
    return pd.DataFrame(out)


def overlap95(a_q: tuple[float, float], b_q: tuple[float, float]) -> bool:
    return not (a_q[1] < b_q[0] or b_q[1] < a_q[0])


def main() -> None:
    if not OLD_LOG.is_file() or not NEW_LOG.is_file():
        sys.exit(f"Missing one of:\n  {OLD_LOG}\n  {NEW_LOG}")

    df_old = load_trace(OLD_LOG)
    df_new = load_trace(NEW_LOG)

    last_old = last_progress(OLD_RUN)
    last_new = last_progress(NEW_RUN)

    print("=" * 70)
    print("Tier-3 OLD vs NEW: walltime + posterior comparison")
    print("=" * 70)
    print(f"OLD: last sample = {last_old[0]:,}, last reported {last_old[1]}/Msample"
          f" (≈ {hms_to_hours(last_old[1]):.2f} h/Msample)")
    print(f"NEW: last sample = {last_new[0]:,}, last reported {last_new[1]}/Msample"
          f" (≈ {hms_to_hours(last_new[1]):.2f} h/Msample)")

    # Wall-time comparison: equal-wall-clock samples-per-time ratio.
    if last_old[0] > 0 and last_new[0] > 0:
        ratio = last_new[0] / last_old[0]
        print(f"NEW reached {last_new[0] / last_old[0]:.3f}× as many samples as OLD "
              f"in equal wall time → speedup ≈ {ratio:.2f}×")

    # Headline parameters to compare.
    headline = [
        "posterior", "likelihood", "prior", "Mascot",
        "Tree.height", "Tree.treeLength",
        "f_migrationRatesSkyline.Deme1_to_Deme2",
        "f_migrationRatesSkyline.Deme2_to_Deme1",
        "f_migrationRatesSkyline.Deme3_to_Deme1",
        "NeScaler.Deme1", "NeScaler.Deme2", "NeScaler.Deme3", "NeScaler.Deme4",
        "wastewater.scaling.Deme1:SimDataset",
        "wastewater.scaling.Deme2:SimDataset",
        "wastewater.scaling.Deme3:SimDataset",
        "caseCounts.dispersion:SimDataset",
    ]
    headline = [h for h in headline if h in df_old.columns and h in df_new.columns]

    s_old = summary(df_old, headline)
    s_new = summary(df_new, headline)
    merged = s_old.merge(s_new, on="param", suffixes=("_old", "_new"))

    print()
    print("Per-parameter summary (post-burn-in 80% of trace, ESS = initial-positive-sequence):")
    print(f"{'param':<45} "
          f"{'old mean':>11} {'old 95%':>22} {'old ESS':>8}  "
          f"{'new mean':>11} {'new 95%':>22} {'new ESS':>8}  overlap?")
    for _, r in merged.iterrows():
        old_q = (r["q025_old"], r["q975_old"])
        new_q = (r["q025_new"], r["q975_new"])
        ov = "YES" if overlap95(old_q, new_q) else "NO "
        print(
            f"{r['param']:<45} "
            f"{r['mean_old']:>11.4g} [{old_q[0]:>9.4g},{old_q[1]:>9.4g}] {r['ess_old']:>8.1f}  "
            f"{r['mean_new']:>11.4g} [{new_q[0]:>9.4g},{new_q[1]:>9.4g}] {r['ess_new']:>8.1f}  {ov}"
        )

    # ESS/hour using equal-wall-time, conservative: assume both ran for same wall
    # clock, take the wall hours/Msample reported.
    print()
    print("Per-parameter ESS / hour of wall clock (post-burn-in):")
    print(f"{'param':<45} {'old ESS/hr':>11} {'new ESS/hr':>11} {'speedup':>9}")
    for _, r in merged.iterrows():
        # Wall hours used for the post-burn-in window:
        h_old = hms_to_hours(last_old[1]) * (last_old[0] * 0.8 / 1e6)
        h_new = hms_to_hours(last_new[1]) * (last_new[0] * 0.8 / 1e6)
        if h_old <= 0 or h_new <= 0:
            continue
        eph_old = r["ess_old"] / h_old
        eph_new = r["ess_new"] / h_new
        ratio = eph_new / eph_old if eph_old > 0 else float("nan")
        print(f"{r['param']:<45} {eph_old:>11.2f} {eph_new:>11.2f} {ratio:>9.2f}×")

    # Trace plot for the headline parameters.
    fig, axes = plt.subplots(4, 2, figsize=(12, 12), sharex=False)
    plot_params = [
        "posterior", "Mascot",
        "Tree.height", "Tree.treeLength",
        "NeScaler.Deme1", "NeScaler.Deme2",
        "f_migrationRatesSkyline.Deme1_to_Deme2",
        "f_migrationRatesSkyline.Deme2_to_Deme1",
    ]
    plot_params = [p for p in plot_params if p in df_old.columns and p in df_new.columns]
    for ax, param in zip(axes.ravel(), plot_params):
        ax.plot(df_old["Sample"], df_old[param], color="tab:red", linewidth=0.6,
                alpha=0.7, label="old")
        ax.plot(df_new["Sample"], df_new[param], color="tab:blue", linewidth=0.6,
                alpha=0.7, label="new")
        ax.set_title(param, fontsize=9)
        ax.legend(loc="best", fontsize=7)
        ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out_png = REPO / "sandbox" / "tier3_traces.png"
    fig.savefig(out_png, dpi=140)
    print(f"\nWrote trace plot: {out_png.relative_to(REPO)}")


if __name__ == "__main__":
    main()
