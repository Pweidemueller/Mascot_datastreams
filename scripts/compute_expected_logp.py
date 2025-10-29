#!/usr/bin/env python3
"""
Compute expected log-likelihood (logP) for case counts using the same logic as
CaseCountLikelihood + NotAKnotSpline in the Java codebase, but independently in Python.

Key points replicated:
- Not-a-knot cubic spline on log-prevalence at knot times (RateShifts times)
- Precompute prevalence on a regular grid (gridRateShifts)
- Lookup prevalence at observation times by snapping to the nearest grid point
- Gamma-Poisson (Negative Binomial) log PMF with mean = prevalence * scaling and dispersion alpha

No external dependencies beyond the Python standard library.
"""

import math
from typing import List

import numpy as np
import scipy as scp

import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute expected log-likelihood (logP) for case counts using the same logic as CaseCountLikelihood + NotAKnotSpline in the Java codebase, but independently in Python."
    )
    parser.add_argument(
        "--rate_shifts",
        type=float,
        nargs="+",
        default=None,
        help="Rate shifts",
    )
    parser.add_argument(
        "--scaling0",
        type=float,
        default=1.0,
        help="Scaling factor for deme 0",
    )
    parser.add_argument(
        "--scaling1",
        type=float,
        default=1.0,
        help="Scaling factor for deme 1",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.5,
        help="Dispersion parameter",
    )
    return parser.parse_args()


def build_natural_cubic_spline(knots_times: List[float], knots_values: List[float]):
    """Create a natural cubic spline S(t) using SciPy."""
    if len(knots_times) != len(knots_values):
        raise ValueError("knots_times and knots_values must have same length")
    if any(knots_times[i] >= knots_times[i + 1] for i in range(len(knots_times) - 1)):
        raise ValueError("knots_times must be strictly increasing")
    return scp.interpolate.make_interp_spline(
        np.asarray(knots_times, dtype=float),
        np.asarray(knots_values, dtype=float),
        k=3,
        bc_type="natural",
    )


def nearest_grid_value(
    grid_times: List[float], grid_values: List[float], t: float
) -> float:
    """Mimic NotAKnotSpline.getPrevalenceAtGridPoint: choose nearest grid point (left on tie)."""
    if t <= grid_times[0]:
        return grid_values[0]
    if t >= grid_times[-1]:
        return grid_values[-1]
    # Find left index such that grid_times[i] <= t < grid_times[i+1]
    lo, hi = 0, len(grid_times) - 2
    while lo <= hi:
        mid = (lo + hi) // 2
        if grid_times[mid] <= t < grid_times[mid + 1]:
            i = mid
            break
        elif t < grid_times[mid]:
            hi = mid - 1
        else:
            lo = mid + 1
    else:
        i = len(grid_times) - 2
    # Choose closer; left on tie
    if abs(grid_times[i] - t) <= abs(grid_times[i + 1] - t):
        return grid_values[i]
    return grid_values[i + 1]


def gamma_poisson_logpmf(k: float, mean: float, alpha: float) -> float:
    """Log PMF of Negative Binomial with mean and dispersion alpha (as in Java).

    r = 1/alpha, p = r / (r + mean)
    log PMF = ln Γ(r + x) - ln Γ(r) - ln Γ(x+1) + r ln p + x ln(1 - p)
    """
    if mean <= 0.0 or alpha <= 0.0 or k < 0.0:
        return float("-inf")
    x = int(round(k))
    r = 1.0 / alpha
    p = r / (r + mean)
    # Clamp as Java code does to avoid log(0)
    p = min(1.0 - 1e-16, max(1e-16, p))
    return (
        math.lgamma(r + x)
        - math.lgamma(r)
        - math.lgamma(x + 1.0)
        + r * math.log(p)
        + x * math.log(1.0 - p)
    )


def main():
    args = parse_args()
    # Inputs mirrored from the updated Java test
    # Observation times and counts per deme
    times0 = [0.0, 0.1, 0.52, 0.98, 1.47, 2.0]
    counts0 = [5.0, 7.0, 6.0, 8.0, 10.0, 0.0]
    times1 = [0.0, 0.15, 0.56, 0.98, 1.78, 1.98]
    counts1 = [3.0, 4.0, 5.0, 6.0, 8.0, 2.0]

    # Knot times (RateShifts value string in test)
    if args.rate_shifts is None:
        rate_shifts = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
    else:
        rate_shifts = [float(i) for i in args.rate_shifts]

    # Grid times (gridRateShifts in test)
    grid_times = [i / 10.0 for i in range(0, 21)]  # 0.0 .. 2.0 step 0.1

    # log-infected control points per knot (already in log space)
    logI_deme0 = [0.0, 1.0, 2.0, 5.0, 10.0, 2.5, -2.0, 0.0, 1.0, -1.0, 0.0]
    logI_deme1 = [0.0, 1.0, 3.0, 6.0, 8.0, 3.0, -2.0, 0.0, -0.5, -1.0, 0.0]

    # Distribution parameters
    alpha = args.alpha
    scaling0 = args.scaling0
    scaling1 = args.scaling1

    # Build natural cubic splines on logI
    spline0 = build_natural_cubic_spline(rate_shifts, logI_deme0)
    spline1 = build_natural_cubic_spline(rate_shifts, logI_deme1)

    # Precompute prevalence on grid: I = exp(spline(logI))
    grid_I0 = [math.exp(float(spline0(t))) for t in grid_times]
    grid_I1 = [math.exp(float(spline1(t))) for t in grid_times]

    # Accumulate log-likelihood
    logP = 0.0
    for t, k in zip(times0, counts0):
        mean = nearest_grid_value(grid_times, grid_I0, t) * scaling0
        logP += gamma_poisson_logpmf(k, mean, alpha)
    for t, k in zip(times1, counts1):
        mean = nearest_grid_value(grid_times, grid_I1, t) * scaling1
        logP += gamma_poisson_logpmf(k, mean, alpha)

    print(f"logP = {logP:.12f}")


if __name__ == "__main__":
    main()
