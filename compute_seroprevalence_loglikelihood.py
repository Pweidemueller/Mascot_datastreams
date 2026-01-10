#!/usr/bin/env python3
"""
Compute log-likelihood for seroprevalence observations using natural spline interpolation.

This script replicates the logic from SeroprevalenceLikelihood.java:
- Creates a natural spline from log-prevalence knot values
- Pre-calculates prevalence and transmission rates at grid points
- For each observation, computes cumulative incidence by integrating incidence rate
- Computes log-likelihood using Binomial distribution

Model:
- Incidence rate = β(t) * I(t) where β is transmission rate and I is prevalence
- Cumulative incidence = integral of incidence rate from earliestTime to observation time
- Seroprevalence probability p = scaling * cumulativeIncidence / populationSize
- Observations follow Binomial(n, p) where n is number tested and x is number seropositive
"""

import numpy as np
import scipy.interpolate
from scipy.special import gammaln
from typing import Optional, Union, List


def binomial_logpmf(x: int, n: int, p: float) -> float:
    """
    Compute log PMF for Binomial distribution: P(X = x) = C(n, x) * p^x * (1-p)^(n-x)

    Parameters
    ----------
    x : int
        Number of successes (seropositive)
    n : int
        Number of trials (people tested)
    p : float
        Success probability (must be in (0, 1))

    Returns
    -------
    float
        Log probability mass function value
    """
    if n < 0:
        return -np.inf
    if x < 0 or x > n:
        return -np.inf

    # Clamp p for numerical stability (matching Java PROB_EPS = 1e-16)
    PROB_EPS = 1e-16
    p = min(1.0 - PROB_EPS, max(PROB_EPS, p))

    # log C(n, x) = log(n!) - log(x!) - log((n-x)!)
    # Using log-gamma: log(n!) = gammaln(n+1)
    log_comb = gammaln(n + 1.0) - gammaln(x + 1.0) - gammaln(n - x + 1.0)

    # log(p^x * (1-p)^(n-x)) = x * log(p) + (n-x) * log(1-p)
    log_p = np.log(p)
    log_1mp = np.log(1.0 - p)

    return log_comb + x * log_p + (n - x) * log_1mp


def compute_cumulative_incidence(
    knots_times: np.ndarray,
    knots_log_prevalence: np.ndarray,
    grid_times: np.ndarray,
    earliest_time: float,
    observation_time: float,
    uninfectious_rate: float,
) -> float:
    """
    Compute cumulative incidence by integrating incidence rate β(t) * I(t) from earliest_time to observation_time.

    This matches the Java getCumulativeIncidence() behavior:
    - Times are in backward time (larger = further in past)
    - Integration uses trapezoidal rule over grid segments
    - Incidence rate = transmission_rate * prevalence
    - Transmission rate β = γ - d(logI)/dτ where γ is uninfectious rate and τ is backward time

    Parameters
    ----------
    knots_times : np.ndarray
        Time points for spline knots (e.g., rate shift times)
    knots_log_prevalence : np.ndarray
        Log-prevalence values at knot points
    grid_times : np.ndarray
        Grid point times where prevalence and transmission rates are pre-calculated
    earliest_time : float
        Start time for integration (further in past, larger value in backward time)
    observation_time : float
        End time for integration (closer to present, smaller value in backward time)
    uninfectious_rate : float
        Rate at which individuals become uninfectious (γ)

    Returns
    -------
    float
        Cumulative incidence (integral of incidence rate)
    """
    # Times are in backward time: if earliest_time <= observation_time, return 0
    if earliest_time <= observation_time:
        return 0.0

    if len(grid_times) <= 1:
        return 0.0

    grid_min = grid_times[0]
    grid_max = grid_times[-1]

    # Clamp to grid range in backward time
    from_time = min(earliest_time, grid_max)
    to_time = max(observation_time, grid_min)

    if from_time <= to_time:
        return 0.0

    # Create natural spline interpolation
    spline = scipy.interpolate.make_interp_spline(
        knots_times, knots_log_prevalence, k=3, bc_type="natural"
    )

    # Pre-calculate prevalence and transmission rates at grid points
    grid_prevalence = np.zeros(len(grid_times))
    grid_transmission_rate = np.zeros(len(grid_times))

    t_min = knots_times[0]
    t_max = knots_times[-1]

    for i, grid_t in enumerate(grid_times):
        # Clamp to boundary values if outside spline range
        if grid_t < t_min:
            logI = spline(t_min)
        elif grid_t > t_max:
            logI = spline(t_max)
        else:
            logI = spline(grid_t)

        grid_prevalence[i] = np.exp(logI)

        # Compute transmission rate: β = γ - d(logI)/dτ
        # where d(logI)/dτ is the derivative of the spline in backward time
        if grid_t < t_min:
            dlogI_dtau = spline.derivative()(t_min)
        elif grid_t > t_max:
            dlogI_dtau = spline.derivative()(t_max)
        else:
            dlogI_dtau = spline.derivative()(grid_t)

        grid_transmission_rate[i] = uninfectious_rate - dlogI_dtau

    # Find boundary indices in backward time
    # getRightGridIndex: smallest i such that grid_times[i] >= from_time
    # getLeftGridIndex: largest i such that grid_times[i] <= to_time

    # Find maxIdx (right grid index for from_time) - matching Java getRightGridIndex
    n_grid = len(grid_times)
    if n_grid == 0:
        return 0.0

    first = grid_times[0]
    last = grid_times[-1]

    if from_time <= first:
        max_idx = 0
    elif from_time >= last:
        max_idx = n_grid - 1
    else:
        # Binary search for smallest i such that grid_times[i] >= from_time
        left = 0
        right = n_grid - 1
        while left + 1 < right:
            mid = (left + right) // 2
            tm = grid_times[mid]
            if tm >= from_time:
                right = mid
            else:
                left = mid
        max_idx = right

    # Find minIdx (left grid index for to_time) - matching Java getLeftGridIndex
    if to_time <= first:
        min_idx = 0
    elif to_time >= last:
        min_idx = n_grid - 1
    else:
        # Binary search for largest i such that grid_times[i] <= to_time
        left = 0
        right = n_grid - 1
        while left + 1 < right:
            mid = (left + right) // 2
            tm = grid_times[mid]
            if tm <= to_time:
                left = mid
            else:
                right = mid
        min_idx = left

    # Integrate using trapezoids over gridpoint segments (backward time: from max_idx down to min_idx)
    sum_incidence = 0.0

    for i in range(max_idx, min_idx, -1):
        t1 = grid_times[i]
        t2 = grid_times[i - 1]

        # Incidence rate = transmission_rate * prevalence
        incidence1 = grid_transmission_rate[i] * grid_prevalence[i]
        incidence2 = grid_transmission_rate[i - 1] * grid_prevalence[i - 1]

        # Trapezoidal rule: (t1 - t2) * (incidence1 + incidence2) * 0.5
        # Note: t1 > t2 in backward time
        sum_incidence += (t1 - t2) * (incidence1 + incidence2) * 0.5

    return sum_incidence


def compute_log_likelihood(
    knots_times: np.ndarray,
    knots_log_prevalence: np.ndarray,
    grid_times: np.ndarray,
    observation_times: np.ndarray,
    people_tested: np.ndarray,
    people_seropositive: np.ndarray,
    population_size: float,
    uninfectious_rate: float,
    scaling: float = 1.0,
    earliest_time: Optional[float] = None,
) -> float:
    """
    Compute log-likelihood for seroprevalence observations.

    Parameters
    ----------
    knots_times : np.ndarray
        Time points for spline knots (e.g., rate shift times)
    knots_log_prevalence : np.ndarray
        Log-prevalence values at knot points
    grid_times : np.ndarray
        Grid point times where prevalence and transmission rates are pre-calculated
    observation_times : np.ndarray
        Times at which seroprevalence was observed
    people_tested : np.ndarray
        Number of people tested at each observation time (must be integers)
    people_seropositive : np.ndarray
        Number of people seropositive at each observation time (must be integers)
    population_size : float
        Population size
    uninfectious_rate : float
        Rate at which individuals become uninfectious (γ)
    scaling : float, optional
        Scaling factor applied to cumulative incidence (default: 1.0)
    earliest_time : float, optional
        Start time for cumulative incidence calculation (default: latest grid time)

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
    if len(observation_times) != len(people_tested) or len(observation_times) != len(
        people_seropositive
    ):
        raise ValueError(
            "observation_times, people_tested, and people_seropositive must have the same length"
        )
    if scaling <= 0:
        raise ValueError("scaling must be > 0")
    if population_size <= 0:
        raise ValueError("population_size must be > 0")
    if uninfectious_rate <= 0:
        raise ValueError("uninfectious_rate must be > 0")

    # Default earliest time is latest grid time (farthest in past)
    if earliest_time is None:
        earliest_time = grid_times[-1]

    # Compute log-likelihood for each observation
    log_likelihood = 0.0

    for t, n, x in zip(observation_times, people_tested, people_seropositive):
        # Validate observation bounds
        n_int = int(np.round(n))
        x_int = int(np.round(x))

        if n_int < 0 or x_int < 0 or x_int > n_int:
            return -np.inf

        # Compute cumulative incidence
        cumulative_incidence = compute_cumulative_incidence(
            knots_times=knots_times,
            knots_log_prevalence=knots_log_prevalence,
            grid_times=grid_times,
            earliest_time=earliest_time,
            observation_time=t,
            uninfectious_rate=uninfectious_rate,
        )

        # Compute seroprevalence probability
        p = scaling * cumulative_incidence / population_size

        # Validate and clamp p
        if np.isnan(p):
            return -np.inf

        # Clamp to (0, 1) for numerical stability (matching Java PROB_EPS = 1e-16)
        PROB_EPS = 1e-16
        p = min(1.0 - PROB_EPS, max(PROB_EPS, p))

        # Compute binomial log-likelihood
        log_pdf = binomial_logpmf(x_int, n_int, p)

        if np.isnan(log_pdf) or np.isinf(log_pdf):
            return -np.inf

        log_likelihood += log_pdf

    return log_likelihood


def main():
    """
    Compute log-likelihoods for test case in SeroprevalenceLikelihoodTest.java
    """

    # Common spline setup
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

    logI = np.array([0.0, 1.0, 2.0, 5.0, 10.0, 2.5, -2.0, 0.0, 1.0, -1.0, 0.0])

    # Test: constantPrevalence_integral_drivesLogLikelihood
    print("Test constantPrevalence_integral_drivesLogLikelihood")

    # One observation at time t = 1.25
    observation_times = np.array([1.25])
    people_tested = np.array([100])
    people_seropositive = np.array([3])

    uninfectious_rate = 1.0  # Transmission rate = 1 for constant prevalence
    population_size = 10000.0
    scaling = 0.5

    # For constant prevalence c and transmission rate = 1, cumulative incidence = c * t
    # Expected: cumulative incidence = 10 * 1.25 = 12.5
    # p = 0.5 * 12.5 / 10000 = 0.00625

    logP = compute_log_likelihood(
        knots_times=knots_times,
        knots_log_prevalence=logI,
        grid_times=grid_times,
        observation_times=observation_times,
        people_tested=people_tested,
        people_seropositive=people_seropositive,
        population_size=population_size,
        uninfectious_rate=uninfectious_rate,
        scaling=scaling,
    )

    print(f"logP = {logP}")


if __name__ == "__main__":
    main()
