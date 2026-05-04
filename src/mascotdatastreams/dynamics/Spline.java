package mascotdatastreams.dynamics;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.RealParameter;
import mascot.dynamics.RateShifts;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;


/**
 * @author Nicola F. Mueller & Paula Weidemueller
 */
@Description("Population function with values at certain time points that are interpolated in between. Parameter has to be in log space. This class uses a natural cubic spline interpolation.")
public class Spline extends CalculationNode {

    final public Input<RealParameter> InfectedInput = new Input<>("logInfected",
            "Prevalence over time in log space", Input.Validate.REQUIRED);
    final public Input<RateShifts> rateShiftsInput = new Input<>("rateShifts",
            "Knots in log prevalence space for spline interpolation, should be given in absolute time units (backward from present)", Input.Validate.REQUIRED);
    final public Input<RateShifts> gridRateShiftsInput = new Input<>("gridRateShifts",
            "Rate shifts to use as grid points for spline evaluation, should be given in absolute time units (backward from present)", Input.Validate.REQUIRED);
    final public Input<RealParameter> uninfectiousRateInput = new Input<>("uninfectiousRate",
            "Rate at which individuals become uninfectious", Input.Validate.REQUIRED);
    final public Input<Boolean> clipTransRateInput = new Input<>("clipTransRate",
            "If true, clip transmission rate to minimum value TR_MIN; if false, use raw rate (default true)",
            true);

    RealParameter infected;
    RateShifts rateShifts;
    RateShifts gridRateShifts;
    RealParameter uninfectiousRate;

    double[] transmissionRate;

    double[] logI;

    PolynomialSplineFunction splineFunction;
    UnivariateFunction splineDerivative;
    PolynomialSplineFunction splineFunction_stored;
    UnivariateFunction splineDerivative_stored;

    double[] time;

    boolean ratesKnows=false;
    boolean isValid = true;
    private boolean clipTransRate;
    private static final double TR_MIN = 1e-1;

    @Override
    public void initAndValidate() {
        infected = InfectedInput.get();
        rateShifts = rateShiftsInput.get();
        gridRateShifts = gridRateShiftsInput.get();
        clipTransRate = clipTransRateInput.get();
        // infected dimension should match rateShifts dimension (rateShifts already includes time 0)
        infected.setDimension(rateShifts.getDimension());
        uninfectiousRate = uninfectiousRateInput.get();
        recalculateRates();
    }

    // Precomputes log-prevalence and transmission rate on the evaluation grid.
    private void recalculateRates() {
        buildSpline();
        // Cache derivative once per recompute: PolynomialSplineFunction#derivative()
        // allocates a new spline, so we don't want to call it per query.
        splineDerivative = splineFunction.derivative();

        int n = gridRateShifts.getDimension();
        time = new double[n];
        logI = new double[n];
        transmissionRate = new double[n];

        isValid = true;

        double firstKnot = rateShifts.getValue(0);
        double lastKnot = rateShifts.getValue(rateShifts.getDimension() - 1);
        double gamma = uninfectiousRate.getValue();

        for (int i = 0; i < n; i++) {
            time[i] = gridRateShifts.getValue(i);

            // Clamp evaluation time to the spline's knot domain.
            // TODO: check if we want to throw an error or take value at first or last knot
            double tEval = time[i];
            if (tEval < firstKnot) tEval = firstKnot;
            else if (tEval > lastKnot) tEval = lastKnot;

            logI[i] = splineFunction.value(tEval);

            // Forward time: d(log I)/dt_fwd = β - γ. Backward time τ flips the sign,
            // so β = γ - d(log I)/dτ.
            double dLogI_dBackwardTime = splineDerivative.value(tEval);
            transmissionRate[i] = gamma - dLogI_dBackwardTime;
            if (clipTransRate) {
                transmissionRate[i] = Math.max(transmissionRate[i], TR_MIN);
            }
        }

        ratesKnows = true;
    }

    public boolean update() {
        if (!ratesKnows) {
            recalculateRates();
        }
        return isValid;
    }

    @Override
    public boolean requiresRecalculation() {
        ratesKnows = false;
        return true;
    }

    @Override
    public void store() {
        splineFunction_stored = splineFunction;
        splineDerivative_stored = splineDerivative;
        super.store();
    }

    @Override
    public void restore() {
        ratesKnows=false;
        splineFunction = splineFunction_stored;
        splineDerivative = splineDerivative_stored;
        super.restore();
    }


    /**
     * Largest index k with time[k] <= t. Assumes time[0] < t < time[last] (callers
     * handle out-of-range clamping).
     */
    private int findLeftIndex(double t) {
        int left = 0, right = time.length - 1;
        while (left < right - 1) {
            int mid = (left + right) / 2;
            if (time[mid] <= t) left = mid; else right = mid;
        }
        return left;
    }

    /**
     * Linear interpolation of a per-grid-point value at time t. Out-of-range
     * t is clamped to the boundary value.
     */
    private double interpolate(double[] values, double t) {
        if (t <= time[0]) return values[0];
        if (t >= time[time.length - 1]) return values[values.length - 1];
        int left = findLeftIndex(t);
        double w = (t - time[left]) / (time[left + 1] - time[left]);
        return values[left] + w * (values[left + 1] - values[left]);
    }

    /**
     * Log-prevalence at time t, linearly interpolated between grid points.
     *
     * @param t time (backward from present)
     */
    public double getLogPrevalence(double t) {
        if (!ratesKnows) recalculateRates();
        return interpolate(logI, t);
    }

    /**
     * Prevalence at time t. Interpolation happens in log space (matching the
     * spline's native parameterisation), then exponentiated.
     *
     * @param t time (backward from present)
     */
    public double getPrevalence(double t) {
        return Math.exp(getLogPrevalence(t));
    }

    /**
     * Transmission rate β(t) = γ - d(logI)/dτ, linearly interpolated between
     * grid points (with optional clipping to TR_MIN applied at grid points).
     *
     * @param t time (backward from present)
     */
    public double getTransmissionRate(double t) {
        if (!ratesKnows) recalculateRates();
        return interpolate(transmissionRate, t);
    }

    /**
     * Exact spline derivative d(log I)/dτ at time t (no grid interpolation).
     * Returns 0 outside the knot range.
     *
     * @param t time (backward from present)
     */
    public double getLogPrevalenceDerivative(double t) {
        if (!ratesKnows) recalculateRates();
        // TODO: check if we want to throw an error or return 0.0
        if (t < time[0] || t > time[time.length - 1]) {
            return 0.0;
        }
        return splineDerivative.value(t);
    }

    /**
     * Get the number of grid points in the spline.
     *
     * @return number of grid points
     */
    public int getGridPointCount() {
        return gridRateShifts.getDimension();
    }

    /**
     * Get the start time of the grid.
     *
     * @return start time
     */
    public double getGridStart() {
        if (gridRateShifts != null && gridRateShifts.getDimension() > 0) {
            return gridRateShifts.getValue(0);
        }
        return 0.0;
    }

    /**
     * Get the end time of the grid.
     *
     * @return end time
     */
    public double getGridEnd() {
        if (gridRateShifts != null && gridRateShifts.getDimension() > 0) {
            return gridRateShifts.getValue(gridRateShifts.getDimension() - 1); // last grid point
        }
        return 0.0;
    }

    public double getGridPointTime(int i) {
        return gridRateShifts.getValue(i);
    }

    /**
     * Returns the largest index i such that grid time[i] <= t.
     * Clamps to [0, getGridPointCount()-1].
     */
    public int getLeftGridIndex(double t) {
        int n = getGridPointCount();
        if (n == 0) return 0;
        double first = getGridPointTime(0);
        double last = getGridPointTime(n - 1);
        if (t <= first) return 0;
        if (t >= last) return n - 1;
        int left = 0;
        int right = n - 1;
        while (left + 1 < right) {
            int mid = (left + right) / 2;
            double tm = getGridPointTime(mid);
            if (tm <= t) {
                left = mid;
            } else {
                right = mid;
            }
        }
        return left;
    }

    /**
     * Returns the smallest index i such that grid time[i] >= t.
     * Clamps to [0, getGridPointCount()-1].
     */
    public int getRightGridIndex(double t) {
        int n = getGridPointCount();
        if (n == 0) return 0;
        double first = getGridPointTime(0);
        double last = getGridPointTime(n - 1);
        if (t <= first) return 0;
        if (t >= last) return n - 1;
        int left = 0;
        int right = n - 1;
        while (left + 1 < right) {
            int mid = (left + right) / 2;
            double tm = getGridPointTime(mid);
            if (tm >= t) {
                right = mid;
            } else {
                left = mid;
            }
        }
        return right;
    }

    /**
     * Get the uninfectious rate parameter.
     *
     * @return uninfectious rate
     */
    public RealParameter getUninfectiousRate() {
        return uninfectiousRate;
    }

    /**
     * Builds a natural cubic spline interpolation using Apache Commons Math.
     * Boundary conditions are natural (second derivative is 0 at the first and last knot)
     */
    private void buildSpline() {
        // rateShifts already includes time 0 as the first value
        // infected[i] corresponds to rateShifts[i] for i = 0, 1, ..., rateShifts.getDimension()-1

        int n = rateShifts.getDimension();

        // Prepare arrays for spline interpolation
        double[] knotTimes = new double[n];
        double[] knotValues = new double[n];

        for (int i = 0; i < n; i++) {
            knotTimes[i] = rateShifts.getValue(i);
            knotValues[i] = infected.getArrayValue(i);
        }

        // Create the spline using Apache Commons Math SplineInterpolator
        // This creates a natural cubic spline (natural boundary conditions)
        SplineInterpolator interpolator = new SplineInterpolator();
        splineFunction = interpolator.interpolate(knotTimes, knotValues);
    }
}
