#!/usr/bin/env python3
"""
Compute log-likelihood for case count observations using natural spline interpolation.

This script replicates the logic from CaseCountLikelihood.java + Spline.java:
- Creates a natural cubic spline from log-prevalence knot values
- Pre-calculates log-prevalence at grid points (matches Java Spline.recalculateRates)
- For each observation time, linearly interpolates logI between bordering grid
  points, then exponentiates to get prevalence (matches Java Spline.getPrevalence
  via getLogPrevalence)
- Computes log-likelihood using a Gamma-Poisson (Negative Binomial) PMF with
  mean = prevalence * scaling and dispersion alpha.
"""

import math
from typing import List

import numpy as np
import scipy.interpolate


def gamma_poisson_logpmf(k: float, mean: float, alpha: float) -> float:
    """Log PMF of Negative Binomial parameterised by mean and dispersion alpha.

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


def compute_log_likelihood(
    knots_times: np.ndarray,
    knots_log_prevalence: np.ndarray,
    grid_times: np.ndarray,
    observation_times: np.ndarray,
    case_counts: np.ndarray,
    alpha: float,
    scaling: float = 1.0,
) -> float:
    """
    Compute log-likelihood for case count observations.

    This matches the Java Spline.getLogPrevalence() / getPrevalence() behaviour:
    - Pre-calculates log-prevalence at grid points
    - For each observation time, linearly interpolates in log space between the
      two bordering grid points (clamping at the boundaries), then exponentiates.

    Parameters
    ----------
    knots_times : np.ndarray
        Time points for spline knots (e.g., rate shift times)
    knots_log_prevalence : np.ndarray
        Log-prevalence values at knot points
    grid_times : np.ndarray
        Grid point times where log-prevalence is pre-calculated
    observation_times : np.ndarray
        Times at which case counts were observed
    case_counts : np.ndarray
        Observed case counts (non-negative)
    alpha : float
        Dispersion parameter of the Gamma-Poisson (Negative Binomial)
    scaling : float, optional
        Scaling factor applied to prevalence-derived mean (default: 1.0)

    Returns
    -------
    float
        Total log-likelihood across all observations
    """
    # Validate inputs
    if len(knots_times) != len(knots_log_prevalence):
        raise ValueError(
            "knots_times and knots_log_prevalence must have the same length"
        )
    if len(observation_times) != len(case_counts):
        raise ValueError(
            "observation_times and case_counts must have the same length"
        )
    if alpha <= 0:
        raise ValueError("alpha must be > 0")
    if scaling <= 0:
        raise ValueError("scaling must be > 0")
    if np.any(case_counts < 0):
        raise ValueError("All case counts must be >= 0")

    # Create natural spline interpolation (cubic spline with natural boundary conditions)
    spline = scipy.interpolate.make_interp_spline(
        knots_times, knots_log_prevalence, k=3, bc_type="natural"
    )

    # Pre-calculate log-prevalence at grid points (matching Java Spline.recalculateRates())
    # Clamp grid times to knot boundaries if needed
    t_min = knots_times[0]
    t_max = knots_times[-1]

    grid_logI = np.zeros(len(grid_times))
    for i, grid_t in enumerate(grid_times):
        # Clamp to boundary values if outside spline range (matching Java behavior)
        if grid_t < t_min:
            grid_logI[i] = spline(t_min)
        elif grid_t > t_max:
            grid_logI[i] = spline(t_max)
        else:
            grid_logI[i] = spline(grid_t)

    # Helper function: linear interpolation in log space (matches Java getLogPrevalence)
    def get_log_prevalence(t: float) -> float:
        """Linearly interpolate logI between bordering grid points; clamp at edges."""
        if t <= grid_times[0]:
            return grid_logI[0]
        if t >= grid_times[-1]:
            return grid_logI[-1]
        # Binary search: largest k with grid_times[k] <= t
        left = 0
        right = len(grid_times) - 1
        while left < right - 1:
            mid = (left + right) // 2
            if grid_times[mid] <= t:
                left = mid
            else:
                right = mid
        w = (t - grid_times[left]) / (grid_times[left + 1] - grid_times[left])
        return grid_logI[left] + w * (grid_logI[left + 1] - grid_logI[left])

    # Compute log-likelihood for each observation
    log_likelihood = 0.0
    for t, k in zip(observation_times, case_counts):
        mean = math.exp(get_log_prevalence(t)) * scaling
        log_pdf = gamma_poisson_logpmf(k, mean, alpha)
        if math.isnan(log_pdf) or math.isinf(log_pdf):
            return -math.inf
        log_likelihood += log_pdf

    return log_likelihood


def main():
    """
    Compute log-likelihoods for all test cases in CaseCountLikelihoodTest.java
    """

    # Common spline setup (same for all tests)
    knots_times = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
    grid_times = np.array(
        [
            0.0,
            0.1,
            0.2,
            0.3,
            0.4,
            0.5,
            0.6,
            0.7,
            0.8,
            0.9,
            1.0,
            1.1,
            1.2,
            1.3,
            1.4,
            1.5,
            1.6,
            1.7,
            1.8,
            1.9,
            2.0,
        ]
    )

    logI_deme0 = np.array([0.0, 1.0, 2.0, 5.0, 10.0, 2.5, -2.0, 0.0, 1.0, -1.0, 0.0])
    logI_deme1 = np.array([0.0, 1.0, 3.0, 6.0, 8.0, 3.0, -2.0, 0.0, -0.5, -1.0, 0.0])

    # Test 1: testSplinePrevalenceTwoDemes
    print("Test testSplinePrevalenceTwoDemes")
    times0 = np.array([0.0, 0.1, 0.52, 0.98, 1.47, 2.0])
    counts0 = np.array([5.0, 7.0, 6.0, 8.0, 10.0, 0.0])
    times1 = np.array([0.0, 0.15, 0.56, 0.98, 1.78, 1.98])
    counts1 = np.array([3.0, 4.0, 5.0, 6.0, 8.0, 2.0])

    alpha = 0.5
    scaling = 1.0

    logP0 = compute_log_likelihood(
        knots_times=knots_times,
        knots_log_prevalence=logI_deme0,
        grid_times=grid_times,
        observation_times=times0,
        case_counts=counts0,
        alpha=alpha,
        scaling=scaling,
    )

    logP1 = compute_log_likelihood(
        knots_times=knots_times,
        knots_log_prevalence=logI_deme1,
        grid_times=grid_times,
        observation_times=times1,
        case_counts=counts1,
        alpha=alpha,
        scaling=scaling,
    )

    logP = logP0 + logP1
    print(f"logP = {logP}")
    print()

    # Test 2: testSplinePrevalenceTwoDemesScaling
    print("Test testSplinePrevalenceTwoDemesScaling")
    # Same observations as test 1; different alpha and per-deme scaling.
    alpha = 0.1
    scaling0 = 0.1
    scaling1 = 0.05

    logP0 = compute_log_likelihood(
        knots_times=knots_times,
        knots_log_prevalence=logI_deme0,
        grid_times=grid_times,
        observation_times=times0,
        case_counts=counts0,
        alpha=alpha,
        scaling=scaling0,
    )

    logP1 = compute_log_likelihood(
        knots_times=knots_times,
        knots_log_prevalence=logI_deme1,
        grid_times=grid_times,
        observation_times=times1,
        case_counts=counts1,
        alpha=alpha,
        scaling=scaling1,
    )

    logP = logP0 + logP1
    print(f"logP = {logP}")
    print()

    # Test 3: testSplinePrevalenceTwoDemesCaseCountsOutsideTree
    print("Test testSplinePrevalenceTwoDemesCaseCountsOutsideTree")
    times0 = np.array([-0.1, 0.0, 0.1, 0.52, 0.98, 1.47, 2.0])
    counts0 = np.array([2.0, 100.0, 5480.0, 1500.0, 560.0, 34.0, 0.0])
    times1 = np.array([-0.1, 0.0, 0.15, 0.56, 0.98, 1.78, 1.98])
    counts1 = np.array([1.0, 124.0, 178.0, 1000.0, 1487.0, 246.0, 2.0])

    alpha = 0.5
    scaling = 1.0

    logP0 = compute_log_likelihood(
        knots_times=knots_times,
        knots_log_prevalence=logI_deme0,
        grid_times=grid_times,
        observation_times=times0,
        case_counts=counts0,
        alpha=alpha,
        scaling=scaling,
    )

    logP1 = compute_log_likelihood(
        knots_times=knots_times,
        knots_log_prevalence=logI_deme1,
        grid_times=grid_times,
        observation_times=times1,
        case_counts=counts1,
        alpha=alpha,
        scaling=scaling,
    )

    logP = logP0 + logP1
    print(f"logP = {logP}")


if __name__ == "__main__":
    main()
