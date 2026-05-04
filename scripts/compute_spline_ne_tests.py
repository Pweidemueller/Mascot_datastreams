#!/usr/bin/env python3
"""
Compute expected values for SplinePrevalenceToNeTest.java.

Mirrors the Java pipeline:
- Build a natural cubic spline on log-prevalence at knot times (RateShifts).
- Precompute logI[i] and transmissionRate[i] on the evaluation grid:
    * tEval is clamped to the knot range [firstKnot, lastKnot].
    * logI[i]            = spline(tEval)
    * transmissionRate[i] = γ - splineDerivative(tEval),   optionally clipped
                           below at TR_MIN = 1e-1 (Java default for clipTransRate).
- For a query time t, linearly interpolate between bordering grid points (clamped
  at the boundaries) — this is what Spline.getLogPrevalence / getPrevalence /
  getTransmissionRate do in Java.
- Ne(t) = NeScaler * I(t) / (coalescentScale * β(t)).

Three test methods are reproduced:
1. testBasicSplineInterpolation  (regular grid, step 0.1)
2. testIrregularGridSplineInterpolation  (irregular grid)
3. testNeSplineInterpolation  (Ne with and without NeScaler)
"""

import math

import numpy as np
import scipy.interpolate


TR_MIN = 1e-1  # Java Spline.TR_MIN, applied when clipTransRate is true (default)


def build_spline_grid(
    knots_times: np.ndarray,
    knots_log_prevalence: np.ndarray,
    grid_times: np.ndarray,
    uninfectious_rate: float,
    clip_trans_rate: bool = True,
):
    """Reproduce Spline.recalculateRates(): return precomputed (logI, transmissionRate)."""
    spline = scipy.interpolate.make_interp_spline(
        knots_times, knots_log_prevalence, k=3, bc_type="natural"
    )
    spline_derivative = spline.derivative()

    first_knot = knots_times[0]
    last_knot = knots_times[-1]

    n = len(grid_times)
    log_I = np.zeros(n)
    transmission_rate = np.zeros(n)
    for i, t in enumerate(grid_times):
        # Clamp evaluation time to knot domain (Java does this in recalculateRates)
        if t < first_knot:
            t_eval = first_knot
        elif t > last_knot:
            t_eval = last_knot
        else:
            t_eval = t
        log_I[i] = spline(t_eval)
        d_log_I = spline_derivative(t_eval)
        beta = uninfectious_rate - d_log_I
        if clip_trans_rate:
            beta = max(beta, TR_MIN)
        transmission_rate[i] = beta

    return log_I, transmission_rate


def interpolate(grid_times: np.ndarray, values: np.ndarray, t: float) -> float:
    """Linear interpolation at time t with boundary clamping (matches Java Spline.interpolate)."""
    if t <= grid_times[0]:
        return float(values[0])
    if t >= grid_times[-1]:
        return float(values[-1])
    # Binary search: largest k with grid_times[k] <= t
    left, right = 0, len(grid_times) - 1
    while left < right - 1:
        mid = (left + right) // 2
        if grid_times[mid] <= t:
            left = mid
        else:
            right = mid
    w = (t - grid_times[left]) / (grid_times[left + 1] - grid_times[left])
    return float(values[left] + w * (values[left + 1] - values[left]))


def get_log_prevalence(grid_times, log_I, t):
    return interpolate(grid_times, log_I, t)


def get_prevalence(grid_times, log_I, t):
    return math.exp(get_log_prevalence(grid_times, log_I, t))


def get_transmission_rate(grid_times, transmission_rate, t):
    return interpolate(grid_times, transmission_rate, t)


def get_ne_time(
    grid_times, log_I, transmission_rate, t, coalescent_scale, ne_scaler=1.0
):
    """SplinePrevalenceToNe.getNeTime: Ne = NeScaler * I / (c * β)."""
    I_t = get_prevalence(grid_times, log_I, t)
    beta = get_transmission_rate(grid_times, transmission_rate, t)
    return ne_scaler * I_t / (coalescent_scale * beta)


def print_assertion(expected, call_str, label, tolerance="tolerance"):
    """Print a Java assertEquals line with the computed expected value."""
    print(f'assertEquals({expected!r}, {call_str}, {tolerance}, "{label}");')


def test_basic_spline_interpolation():
    print("=" * 70)
    print("Test testBasicSplineInterpolation")
    print("=" * 70)

    knots_times = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
    grid_times = np.arange(0.0, 2.0 + 1e-9, 0.1)
    log_I_knots = np.array(
        [0.0, 1.0, 2.0, 5.0, 10.0, 2.5, -2.0, 0.0, 1.0, -1.0, 0.0]
    )
    gamma = 75.0

    log_I, _ = build_spline_grid(knots_times, log_I_knots, grid_times, gamma)

    # Knot points: at knots, the spline interpolates exactly, so log(prevalence) ==
    # the knot value provided. (Both 0.6 and 0.8 are also grid points — exact.)
    for t, expected_log in [
        (0.0, 0.0),
        (0.6, 5.0),
        (0.8, 10.0),
        (1.0, 2.5),
        (1.8, -1.0),
        (2.0, 0.0),
    ]:
        got = math.log(get_prevalence(grid_times, log_I, t))
        print_assertion(
            got, f"Math.log(dyn.getPrevalenceTime({t}))", f"Prevalence at t={t}"
        )

    # Off-knot grid points (these too are exact: returns logI[i] directly)
    for t in [0.1, 1.5]:
        got = get_prevalence(grid_times, log_I, t)
        print_assertion(
            got, f"dyn.getPrevalenceTime({t})", f"Prevalence at t={t}"
        )

    # Off-grid query times (linearly interpolated between bordering grid points)
    for t in [0.01, 0.05, 0.09, 1.91, 1.95, 1.98]:
        got = get_prevalence(grid_times, log_I, t)
        print_assertion(
            got, f"dyn.getPrevalenceTime({t})", f"Prevalence at t={t}"
        )

    # Outside knot range — clamped to nearest grid boundary
    for t in [-0.1, 2.05]:
        got = get_prevalence(grid_times, log_I, t)
        print_assertion(
            got, f"dyn.getPrevalenceTime({t})", f"Prevalence at t={t}"
        )


def test_irregular_grid_spline_interpolation():
    print()
    print("=" * 70)
    print("Test testIrregularGridSplineInterpolation")
    print("=" * 70)

    knots_times = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
    grid_times = np.array([0.0, 0.24, 0.76, 1.68, 2.0])
    log_I_knots = np.array(
        [0.0, 1.0, 2.0, 5.0, 10.0, 2.5, -2.0, 0.0, 1.0, -1.0, 0.0]
    )
    gamma = 75.0

    log_I, _ = build_spline_grid(knots_times, log_I_knots, grid_times, gamma)

    # At grid points: prevalence == exp(spline(t))
    for t in [0.24, 0.76, 1.68]:
        got = get_prevalence(grid_times, log_I, t)
        print_assertion(
            got, f"dyn.getPrevalenceTime({t})", f"Prevalence at t={t}"
        )

    # Off-grid times — linear interpolation in log space between bordering grid points
    for t in [0.01, 0.05, 0.09, 0.15, 0.49, 0.27]:
        got = get_prevalence(grid_times, log_I, t)
        print_assertion(
            got, f"dyn.getPrevalenceTime({t})", f"Prevalence at t={t}"
        )


def test_ne_spline_interpolation():
    print()
    print("=" * 70)
    print("Test testNeSplineInterpolation")
    print("=" * 70)

    knots_times = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
    grid_times = np.arange(0.0, 2.0 + 1e-9, 0.1)
    log_I_knots = np.array(
        [0.0, 1.0, 2.0, 5.0, 10.0, 2.5, -2.0, 0.0, 1.0, -1.0, 0.0]
    )
    gamma = 75.0
    coalescent_scale = 2.0

    log_I, transmission_rate = build_spline_grid(
        knots_times, log_I_knots, grid_times, gamma
    )

    # Without NeScaler (Ne)
    for t in [0.0, 0.1, 1.5, 2.0]:
        got = get_ne_time(
            grid_times, log_I, transmission_rate, t, coalescent_scale
        )
        print_assertion(got, f"dyn.getNeTime({t})", f"Ne at t={t}")

    # Off-grid times — linear interp on logI and transmissionRate, then Ne formula
    for t in [0.48, 1.82]:
        got = get_ne_time(
            grid_times, log_I, transmission_rate, t, coalescent_scale
        )
        print_assertion(got, f"dyn.getNeTime({t})", f"Ne at t={t}")

    # Outside knot intervals — clamped to boundary
    for t in [-0.1, 2.05]:
        got = get_ne_time(
            grid_times, log_I, transmission_rate, t, coalescent_scale
        )
        print_assertion(got, f"dyn.getNeTime({t})", f"Ne at t={t}")

    # With NeScaler = 2.0 (Ne values scale multiplicatively)
    print()
    print("# With NeScaler = 2.0")
    ne_scaler = 2.0
    for t in [0.0, 0.1, 1.5, 2.0]:
        got = get_ne_time(
            grid_times,
            log_I,
            transmission_rate,
            t,
            coalescent_scale,
            ne_scaler=ne_scaler,
        )
        print_assertion(got, f"dynScaled.getNeTime({t})", f"NeScaler at t={t}")


def main():
    test_basic_spline_interpolation()
    test_irregular_grid_spline_interpolation()
    test_ne_spline_interpolation()


if __name__ == "__main__":
    main()
