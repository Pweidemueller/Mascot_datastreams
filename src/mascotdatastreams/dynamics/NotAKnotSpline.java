package mascotdatastreams.dynamics;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.RealParameter;
import mascot.dynamics.RateShifts;
import org.apache.commons.math3.analysis.interpolation.AkimaSplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;


/**
 * @author Nicola F. Mueller & Paula Weidemueller
 */
@Description("Population function with values at certain time points that are interpolated in between. Parameter has to be in log space")
public class NotAKnotSpline extends CalculationNode {

    final public Input<RealParameter> InfectedInput = new Input<>("logInfected",
            "Prevalence over time in log space", Input.Validate.REQUIRED);
    final public Input<RateShifts> rateShiftsInput = new Input<>("rateShifts",
            "Knots in log prevalence space for spline interpolation, should be given in absolute time units (backward from present)", Input.Validate.REQUIRED);
    final public Input<RateShifts> gridRateShiftsInput = new Input<>("gridRateShifts",
            "Rate shifts to use as grid points for spline evaluation, should be given in absolute time units (backward from present)", Input.Validate.REQUIRED);
    final public Input<RealParameter> uninfectiousRateInput = new Input<>("uninfectiousRate",
            "Rate at which individuals become uninfectious", Input.Validate.REQUIRED);

    RealParameter infected;
    RateShifts rateShifts;
    RateShifts gridRateShifts;
    RealParameter uninfectiousRate;

    double[] transmissionRate;
    double[] transmissionRateStored;

    double[] I;
    double[] I_stored;

    PolynomialSplineFunction splineFunction;
    PolynomialSplineFunction splineFunction_stored;

    double[] time;

    boolean ratesKnows=false;
    boolean isValid = true;

    @Override
    public void initAndValidate() {
        infected = InfectedInput.get();
        rateShifts = rateShiftsInput.get();
        gridRateShifts = gridRateShiftsInput.get();
        infected.setDimension(rateShifts.getDimension()+1); // +1 for time 0
        uninfectiousRate = uninfectiousRateInput.get();
        recalculateRates();
    }

    // computes the Ne's at the break points from the growth rates and the transmission rates
    private void recalculateRates() {
        // Build the spline using Apache Commons Math
        buildSpline();

        // use grid rate shifts as grid points
        int n = gridRateShifts.getDimension();
        time = new double[n];
        I = new double[n];
        transmissionRate = new double[n];

        isValid = true;

        // Evaluate spline at all grid points
        for (int i = 0; i < gridRateShifts.getDimension(); i++) {
            time[i] = gridRateShifts.getValue(i);

            // Evaluate spline value (log prevalence)
            double logI = splineFunction.value(time[i]);
            I[i] = Math.exp(logI);

            // Evaluate spline derivative (d(log I)/dt)
            double derivative_logprev = splineFunction.derivative().value(time[i]);

            // Renewal equation: dI/dt = beta*I - gamma*I
            // Since dI/dt = I * d(log I)/dt, we have:
            // I * d(log I)/dt = beta*I - gamma*I
            // d(log I)/dt = beta - gamma
            // beta = gamma + d(log I)/dt
            transmissionRate[i] = uninfectiousRate.getValue() + derivative_logprev;
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
        transmissionRateStored = new double[transmissionRate.length];
        System.arraycopy(transmissionRate, 0, transmissionRateStored, 0, transmissionRate.length);
        I_stored = new double[I.length];
        System.arraycopy(I, 0, I_stored, 0, I.length);
        splineFunction_stored = splineFunction;
        super.store();
    }

    @Override
    public void restore() {
        ratesKnows=false;
        splineFunction = splineFunction_stored;
        super.restore();
    }


    /**
     * Get the log-prevalence value at a specific time using precomputed grid points.
     *
     * @param t time (backward from present)
     * @return log-prevalence at time t
     */
    public double getValueAtGridPoint(double t) {
        if (!ratesKnows) {
            recalculateRates();
        }

        // Find the appropriate grid point
        if (t <= time[0]) {
            return Math.log(I[0]);
        }

        // Binary search for the closest grid point
        int left = 0, right = time.length - 1;
        while (left < right) {
            int mid = (left + right) / 2;
            if (mid + 1 < time.length && time[mid] <= t && time[mid+1] > t) {
                return Math.log(I[mid]);
            }
            if (time[mid] <= t) {
                left = mid+1;
            } else {
                right = mid;
            }
        }

        return Math.log(I[left]);
    }

    /**
     * Get the log-prevalence value at a specific time using precomputed grid points.
     *
     * @param t time (backward from present)
     * @return log-prevalence at time t
     */
    public double getPrevalenceAtGridPoint(double t) {
        if (!ratesKnows) {
            recalculateRates();
        }

        // Find the appropriate grid point
        if (t <= time[0]) {
            return I[0];
        }
        // Binary search for the closest grid point
        int left = 0, right = time.length - 1;
        while (left < right) {
            int mid = (left + right) / 2;
            if (mid + 1 < time.length && time[mid] <= t && time[mid+1] > t) {
                return I[mid];
            }
            if (time[mid] <= t) {
                left = mid+1;
            } else {
                right = mid;
            }
        }
        return I[left];
    }

    /**
     * Get the derivative of log-prevalence at a specific time.
     *
     * @param t time (backward from present)
     * @return derivative of log-prevalence at time t
     */
    public double getDerivativeAtGridPoint(double t) {
        if (!ratesKnows) {
            recalculateRates();
        }

        return splineFunction.derivative().value(t);
    }

    public double getTranssmissionRateAtGridPoint(double t) {
        if (!ratesKnows) {
            recalculateRates();
        }

        // Find the appropriate grid point
        if (t <= time[0]) {
            return transmissionRate[0];
        }

        // Binary search for the closest grid point
        int left = 0, right = time.length - 1;
        while (left < right) {
            int mid = (left + right) / 2;
            if (mid + 1 < time.length && time[mid] <= t && time[mid+1] > t) {
                return transmissionRate[mid];
            }
            if (time[mid] <= t) {
                left = mid+1;
            } else {
                right = mid;
            }
        }

        return transmissionRate[left];
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
        return 0.0;
    }

    /**
     * Get the end time of the grid.
     *
     * @return end time
     */
    public double getGridEnd() {
        if (gridRateShifts != null && gridRateShifts.getDimension() > 0) {
            return gridRateShifts.getValue(gridRateShifts.getDimension() - 1);
        }
        return 0.0;
    }

    public double getGridPointTime(int i) {
        return gridRateShifts.getValue(i);
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
     * Builds the Akima spline interpolation using Apache Commons Math.
     * Akima spline is a C1 differentiable piecewise cubic polynomial that is more stable
     * than cubic splines and does not exhibit overshoot.
     */
    private void buildSpline() {
        int n = rateShifts.getDimension();

        // Prepare arrays for spline interpolation
        double[] knotTimes = new double[n];
        double[] knotValues = new double[n];

        for (int i = 0; i < n; i++) {
            knotTimes[i] = rateShifts.getValue(i);
            knotValues[i] = infected.getArrayValue(i);
        }

        // Create the spline using Apache Commons Math Akima interpolator
        AkimaSplineInterpolator interpolator = new AkimaSplineInterpolator();
        splineFunction = interpolator.interpolate(knotTimes, knotValues);
    }
}
