#!/usr/bin/env python3
"""
Compute log-likelihood for wastewater concentration observations using natural spline interpolation.

This script replicates the logic from WastewaterLikelihood.java:
- Creates a natural spline from log-prevalence knot values
- For each observation, interpolates log-prevalence at observation time
- Computes log-likelihood using LogNormal distribution parameterized by median in real space

Model: log(concentration) ~ Normal(μ, σ²),  where μ = log(α · I(t) / N)
Equivalently: concentration ~ LogNormal(μ, σ²) with median = α · I(t) / N.
I(t) = exp(spline(t)) is the prevalence, N is the population size, α is an optional scaling factor.
"""

import numpy as np
import scipy.interpolate
from scipy.stats import lognorm
from typing import Optional, Union, List


def compute_log_likelihood(
    knots_times: np.ndarray,
    knots_log_prevalence: np.ndarray,
    grid_times: np.ndarray,
    observation_times: np.ndarray,
    concentrations: np.ndarray,
    sd_log: float,
    population_size: float,
    scaling: float = 1.0,
) -> float:
    """
    Compute log-likelihood for wastewater concentration observations.

    This matches the Java Spline.getValueAtGridPoint() behavior:
    - Pre-calculates log-prevalence at grid points
    - For each observation time, returns the closest pre-calculated grid point value

    Parameters
    ----------
    knots_times : np.ndarray
        Time points for spline knots (e.g., rate shift times)
    knots_log_prevalence : np.ndarray
        Log-prevalence values at knot points
    grid_times : np.ndarray
        Grid point times where log-prevalence will be pre-calculated
    observation_times : np.ndarray
        Times at which wastewater concentrations were observed
    concentrations : np.ndarray
        Observed PMV-normalized wastewater concentrations (must be > 0)
    sd_log : float
        Standard deviation on log scale for the LogNormal distribution
    population_size : float
        Population size (used to convert absolute prevalence to per-capita prevalence)
    scaling : float, optional
        Scaling factor applied to prevalence-derived median (default: 1.0)

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
    if len(observation_times) != len(concentrations):
        raise ValueError(
            "observation_times and concentrations must have the same length"
        )
    if sd_log <= 0:
        raise ValueError("sd_log must be > 0")
    if scaling <= 0:
        raise ValueError("scaling must be > 0")
    if population_size <= 0:
        raise ValueError("population_size must be > 0")
    if np.any(concentrations <= 0):
        raise ValueError("All concentrations must be > 0")

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

    # Helper function to find closest grid point (matching Java getValueAtGridPoint)
    def get_value_at_grid_point(t: float) -> float:
        """Find closest pre-calculated grid point value for time t."""
        if t <= grid_times[0]:
            return grid_logI[0]

        # Binary search for the segment containing t
        left = 0
        right = len(grid_times) - 1
        while left < right:
            mid = (left + right) // 2
            if grid_times[mid] <= t and (
                mid == len(grid_times) - 1 or grid_times[mid + 1] > t
            ):
                left = mid
                break
            if grid_times[mid] <= t:
                left = mid + 1
            else:
                right = mid

        # Check if left or right grid point is closer to t
        if left == len(grid_times) - 1:
            return grid_logI[left]
        elif abs(grid_times[left] - t) <= abs(grid_times[left + 1] - t):
            return grid_logI[left]
        else:
            return grid_logI[left + 1]

    # Compute log-likelihood for each observation
    log_likelihood = 0.0

    for t, concentration in zip(observation_times, concentrations):
        # Get log-prevalence at observation time using closest grid point
        logI = get_value_at_grid_point(t)

        # Convert to prevalence (absolute number of infected)
        I = np.exp(logI)

        # Convert to per-capita prevalence (proportion of population infected)
        I_N = I / population_size

        # Apply scaling factor; clip to detection-limit floor (matching Java)
        scaled_median = max(I_N * scaling, 1e-3)

        # Parameterize LogNormal by median in real space.
        # scipy.stats.lognorm(s=σ, scale=exp(μ)): median = exp(μ), so scale = scaled_median.
        log_pdf = lognorm.logpdf(concentration, s=sd_log, scale=scaled_median)

        if np.isnan(log_pdf) or np.isinf(log_pdf):
            return -np.inf

        log_likelihood += log_pdf

    return log_likelihood


def main():
    """
    Compute log-likelihoods for all test cases in WastewaterLikelihoodTest.java
    """

    tree_height = 2.0

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

    # Test 1: testSplinePrevalenceTwoDemes
    print("Test testSplinePrevalenceTwoDemes")
    logI_deme0 = np.array([0.0, 1.0, 2.0, 5.0, 10.0, 2.5, -2.0, 0.0, 1.0, -1.0, 0.0])
    logI_deme1 = np.array([0.0, 1.0, 3.0, 6.0, 8.0, 3.0, -2.0, 0.0, -0.5, -1.0, 0.0])

    times0 = np.array([0.0, 0.1, 0.52, 0.98, 1.47, 2.0])
    concentrations0 = np.array([0.0025, 0.2, 1.5, 35.7, 14.8, 0.8])
    times1 = np.array([0.0, 0.15, 0.56, 0.98, 1.78, 1.98])
    concentrations1 = np.array([0.01, 0.54, 2.3, 4.5, 1.2, 0.9])

    sd_log = 0.5
    scaling = 1.0
    population_size = 10000.0

    logP0 = compute_log_likelihood(
        knots_times=knots_times,
        knots_log_prevalence=logI_deme0,
        grid_times=grid_times,
        observation_times=times0,
        concentrations=concentrations0,
        sd_log=sd_log,
        population_size=population_size,
        scaling=scaling,
    )

    logP1 = compute_log_likelihood(
        knots_times=knots_times,
        knots_log_prevalence=logI_deme1,
        grid_times=grid_times,
        observation_times=times1,
        concentrations=concentrations1,
        sd_log=sd_log,
        population_size=population_size,
        scaling=scaling,
    )

    logP = logP0 + logP1
    print(f"logP = {logP}")
    print()

    # Test 2: testSplinePrevalenceTwoDemesScaling
    print("Test testSplinePrevalenceTwoDemesScaling")
    # Same spline and observations as test 1
    sd_log = 0.1
    population_size = 10000.0
    scaling0 = 0.1
    scaling1 = 0.05

    logP0 = compute_log_likelihood(
        knots_times=knots_times,
        knots_log_prevalence=logI_deme0,
        grid_times=grid_times,
        observation_times=times0,
        concentrations=concentrations0,
        sd_log=sd_log,
        population_size=population_size,
        scaling=scaling0,
    )

    logP1 = compute_log_likelihood(
        knots_times=knots_times,
        knots_log_prevalence=logI_deme1,
        grid_times=grid_times,
        observation_times=times1,
        concentrations=concentrations1,
        sd_log=sd_log,
        population_size=population_size,
        scaling=scaling1,
    )

    logP = logP0 + logP1
    print(f"logP = {logP}")
    print()

    # Test 3: testSplinePrevalenceTwoDemesConcentrationsOutsideTree
    print("Test testSplinePrevalenceTwoDemesConcentrationsOutsideTree")
    # Same spline setup
    times0 = np.array([-0.1, 0.0, 0.1, 0.52, 0.98, 1.47, 2.0])
    concentrations0 = np.array([0.002, 0.0025, 0.2, 1.5, 35.7, 14.8, 0.8])
    times1 = np.array([-0.1, 0.0, 0.15, 0.56, 0.98, 1.78, 1.98])
    concentrations1 = np.array([0.12, 0.01, 0.54, 2.3, 4.5, 1.2, 0.9])

    sd_log = 0.5
    population_size = 10000.0  # Default population size (can be adjusted per test)
    scaling = 1.0

    logP0 = compute_log_likelihood(
        knots_times=knots_times,
        knots_log_prevalence=logI_deme0,
        grid_times=grid_times,
        observation_times=times0,
        concentrations=concentrations0,
        sd_log=sd_log,
        population_size=population_size,
        scaling=scaling,
    )

    logP1 = compute_log_likelihood(
        knots_times=knots_times,
        knots_log_prevalence=logI_deme1,
        grid_times=grid_times,
        observation_times=times1,
        concentrations=concentrations1,
        sd_log=sd_log,
        population_size=population_size,
        scaling=scaling,
    )

    logP = logP0 + logP1
    print(f"logP = {logP}")


if __name__ == "__main__":
    main()
