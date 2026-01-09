package mascotdatastreams.dynamics;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.RealParameter;
import mascot.dynamics.RateShifts;
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

    // TODO revisit this, since ideally we don't need to use clipping but the sampler should reject unreasonable transmissionrate values
    // private static final double TSR_MIN = 0.5;

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
        // infected dimension should match rateShifts dimension (rateShifts already includes time 0)
        infected.setDimension(rateShifts.getDimension());
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

            double logI;
            // TODO: check if we want to throw an error or take value at first or last knot
            if (time[i] < rateShifts.getValue(0)) {
                logI = splineFunction.value(rateShifts.getValue(0)); 
            } else if (time[i] > rateShifts.getValue(rateShifts.getDimension() - 1)) {
                logI = splineFunction.value(rateShifts.getValue(rateShifts.getDimension() - 1)); // logI at last knot
            } else {
                logI = splineFunction.value(time[i]);
            }
            I[i] = Math.exp(logI);

            // Evaluate spline derivative d(log I)/dτ where τ is backward time
            double dLogI_dBackwardTime = splineFunction.derivative().value(time[i]);

            // Epidemiological model in forward time: dI/dt_fwd = (β - γ) * I
            // This means: d(log I)/dt_fwd = β - γ
            // In backward time τ: d(log I)/dτ = -d(log I)/dt_fwd = -(β - γ) = γ - β
            // Rearranging: β = γ - d(log I)/dτ
            // Therefore: transmission_rate = β = γ - d(log I)/dτ_backward
            transmissionRate[i] = uninfectiousRate.getValue() - dLogI_dBackwardTime;
            // TODO: revisit clamp transmission rate to minimum value to prevent division by zero
            // transmissionRate[i] = Math.max(transmissionRate[i], TSR_MIN);
        }

        ratesKnows = true;
    }

    private double findSplineSegmentStartKnotTime(int segmentIndex) {
        return rateShifts.getValue(segmentIndex);
    }

    private double findSplineSegmentEndKnotTime(int segmentIndex) {
        return rateShifts.getValue(segmentIndex+1);
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
        
        // Binary search for the closest sgment in the grid, returns the left grid point of the segment in which t lies
        int left = 0, right = time.length - 1;
        while (left < right) {
            int mid = (left + right ) / 2;
            if (time[mid] <= t && time[mid+1] > t) {
            	left = mid;
            	break;
            }
            if (time[mid] <= t) {
                left = mid+1;
            } else {
                right = mid;
            }
        }
        // check if left or right grid point is closer to t
        if (left == time.length - 1) {
            return Math.log(I[left]);
        }
        else if (Math.abs(time[left] - t) <= Math.abs(time[left+1] - t)) {
        	return Math.log(I[left]);
            
        } 
        return Math.log(I[left+1]);
        
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
        // Binary search for the closest sgment in the grid, returns the left grid point of the segment in which t lies
        int left = 0, right = time.length - 1;
        while (left < right) {
            int mid = (left + right ) / 2;
            if (time[mid] <= t && time[mid+1] > t) {
            	left = mid;
            	break;
            }
            if (time[mid] <= t) {
                left = mid+1;
            } else {
                right = mid;
            }
        }
        // check if left or right grid point is closer to t
        if (left == time.length - 1) {
            return I[left];
        }
        else if (Math.abs(time[left] - t) <= Math.abs(time[left+1] - t)) {
        	return I[left];
            
        } 
        return I[left+1];
    }

    /**
     * Get the derivative of log-prevalence at a specific time.
     *
     * @param t time (backward from present)
     * @return log-prevalence at time t
     */
    public double getDerivative(double t) {
        if (!ratesKnows) {
            recalculateRates();
        }
        // TODO: check if we want to throw an error or return 0.0
        if (t < time[0] || t > time[time.length - 1]) {
            return 0.0;
        } else {
            return splineFunction.derivative().value(t);
        }
    }

    public double getTranssmissionRateAtGridPoint(double t) {
        if (!ratesKnows) {
            recalculateRates();
        }
        // Find the appropriate grid point
        if (t <= time[0]) {
            // Return transmission rate at first point using spline coefficients
            // This is the first knot so dt = 0 
            return transmissionRate[0];
        }
        
         // Binary search for the closest sgment in the grid, returns the left grid point of the segment in which t lies
        int left = 0, right = time.length - 1;
        while (left < right) {
            int mid = (left + right ) / 2;
            if (time[mid] <= t && time[mid+1] > t) {
            	left = mid;
            	break;
            }
            if (time[mid] <= t) {
                left = mid+1;
            } else {
                right = mid;
            }
        }
        // check if left or right grid point is closer to t
        if (left == time.length - 1) {
            return transmissionRate[left];
        }
        else if (Math.abs(time[left] - t) <= Math.abs(time[left+1] - t)) {
        	return transmissionRate[left];
            
        } 
        return transmissionRate[left+1];
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
