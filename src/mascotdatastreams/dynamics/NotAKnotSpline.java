package mascotdatastreams.dynamics;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.RealParameter;
import mascot.dynamics.RateShifts;
import org.apache.commons.math3.linear.*;


/**
 * @author Nicola F. Mueller & Paula Weidemueller
 */
@Description("Populaiton function with values at certain time points that are interpolated in between. Parameter has to be in log space")
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

    double[][] splineCoeffs;
    double[][] splineCoeffs_stored;

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
        notAKnotCubicSpline();
        // use grid rate shifts as grid points
        int n = gridRateShifts.getDimension();
        time = new double[n];
        I = new double[n];
        transmissionRate = new double[n];
        
        isValid = true;

        // First point at time 0
        time[0] = 0.0;
        I[0] = Math.exp(infected.getArrayValue(0));
        double derivative_logprev = splineCoeffs[0][2];
        transmissionRate[0] = uninfectiousRate.getValue() - derivative_logprev; // derivative at t=0 is splineCoeffs[0][2]
        if (transmissionRate[0] < 0) {
            isValid = false;
        }
        
        // Points at grid rate shifts
        for (int i = 1; i < gridRateShifts.getDimension(); i++) {
            time[i] = gridRateShifts.getValue(i);
            // Find which spline segment this grid point falls into
            int segmentIndex = findSplineSegment(gridRateShifts.getValue(i));

            // Bounds check for segmentIndex
            if (segmentIndex >= splineCoeffs.length ) {
                segmentIndex = splineCoeffs.length - 1;
            }
                        
            double timeDiff = gridRateShifts.getValue(i) - findSplineSegmentStartKnotTime(segmentIndex);
            double timeDiff2 = timeDiff * timeDiff;
            double timeDiff3 = timeDiff2 * timeDiff;
            
            // Evaluate spline at this grid point
            I[i] = Math.exp(splineCoeffs[segmentIndex][0]*timeDiff3 + splineCoeffs[segmentIndex][1]*timeDiff2 + 
                             splineCoeffs[segmentIndex][2]*timeDiff + splineCoeffs[segmentIndex][3]);
            
            // Calculate derivative at this grid point
            derivative_logprev = 3 * splineCoeffs[segmentIndex][0] * timeDiff2 + 2 * splineCoeffs[segmentIndex][1] * timeDiff + splineCoeffs[segmentIndex][2];
            transmissionRate[i] = uninfectiousRate.getValue() - derivative_logprev;
            
            if (transmissionRate[i] < 0) {
                isValid = false;
            }
        }
        
        ratesKnows = true;
    }
    
    // Helper method to find which spline segment a given time falls into
    private int findSplineSegment(double t) {
        for (int i = 1; i < rateShifts.getDimension(); i++) {
            if (t < rateShifts.getValue(i)) {
                return i-1;
            }
        }
        // Use last segment if t is beyond all rate shifts
        // splineCoeffs has rateShifts.getDimension() if rateShifts do NOT include time 0.0, otherwise rateShifts.getDimension() - 1 elements (indices 0 to rateShifts.getDimension()-2)
        return rateShifts.getDimension() - 2;
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
        super.store();
    }

    @Override
    public void restore() {
        ratesKnows=false;
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
            int mid = (left + right ) / 2;
            if (time[mid] <= t && time[mid+1] > t) {
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
            int mid = (left + right ) / 2;
            if (time[mid] <= t && time[mid+1] > t) {
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
     * Get the derivative of log-prevalence at a specific time using precomputed grid points.
     * 
     * @param t time (backward from present)
     * @return derivative of log-prevalence at time t
     */
    public double getDerivativeAtGridPoint(double t) {
        if (!ratesKnows) {
            recalculateRates();
        }
        
        // Find the appropriate grid point
        if (t <= time[0]) {
            // Return derivative at first point using spline coefficients
            // This is the first knot so dt = 0 
            return splineCoeffs[0][2];
        }
        
        // Binary search for the closest grid point
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
        
        // Calculate derivative using spline coefficients
        int k = left;
        if (k >= splineCoeffs.length) {
            k = splineCoeffs.length-1;
        }
        double timeDiff = t - findSplineSegmentStartKnotTime(k);
        double timeDiff2 = timeDiff * timeDiff;
        return 3 * splineCoeffs[k][0] * timeDiff2 + 2 * splineCoeffs[k][1] * timeDiff + splineCoeffs[k][2];
    }

    public double getTranssmissionRateAtGridPoint(double t) {
        if (!ratesKnows) {
            recalculateRates();
        }
        // Find the appropriate grid point
        if (t <= time[0]) {
            // Return derivative at first point using spline coefficients
            // This is the first knot so dt = 0 
            return transmissionRate[0];
        }
        
        // Binary search for the closest grid point
        int left = 0, right = time.length - 1;
        while (left < right) {
            int mid = (left + right ) / 2;
            if (time[mid] <= t && time[mid+1] > t) {
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

    /** computes the coefficients for the cubic spline interpolation
     *  de Boor, Carl. A Practical Guide to Splines. Springer-Verlag, New York: 1978
     */
    public void notAKnotCubicSpline() {
        int n = rateShifts.getDimension();

        // Handle edge case: only 2 points - linear interpolation
        if (n == 2) {
            double h = rateShifts.getValue(1) - rateShifts.getValue(0);
            double delta = (infected.getArrayValue(1) - infected.getArrayValue(0)) / h;
            splineCoeffs = new double[1][4];
            splineCoeffs[0][0] = 0.0;  // d = 0 (no cubic term)
            splineCoeffs[0][1] = 0.0;  // c = 0 (no quadratic term)
            splineCoeffs[0][2] = delta; // b = slope
            splineCoeffs[0][3] = infected.getArrayValue(0);  // a = y[0]
            return;
        }

        // Calculate h values (difference between knots)
        double[] h = new double[n - 1];
        for (int i = 0; i < n - 1; i++) {
            h[i] = rateShifts.getValue(i+1) - rateShifts.getValue(i);
        }

        // Calculate the difference in y values
        double[] delta = new double[n - 1];
        for (int i = 0; i < n - 1; i++) {
            delta[i] = (infected.getArrayValue(i + 1) - infected.getArrayValue(i)) / h[i];
        }

        // Create the tridiagonal system
        RealMatrix A = new Array2DRowRealMatrix(n, n);
        RealVector r = new ArrayRealVector(n);

        if (n == 3) {
            // Special case: 3 points - use natural spline boundary conditions
            // Natural spline: second derivative = 0 at endpoints
            // This gives us mu[0] = 0 and mu[n-1] = 0
            
            // First row: natural boundary condition at x[0] (mu[0] = 0)
            A.setEntry(0, 0, 1.0);
            A.setEntry(0, 1, 0.0);
            A.setEntry(0, 2, 0.0);
            r.setEntry(0, 0.0);
            
            // Middle row: continuity condition at x[1]
            A.setEntry(1, 0, h[0]);
            A.setEntry(1, 1, 2 * (h[0] + h[1]));
            A.setEntry(1, 2, h[1]);
            r.setEntry(1, 6 * (delta[1] - delta[0]));
            
            // Last row: natural boundary condition at x[2] (mu[2] = 0)
            A.setEntry(2, 0, 0.0);
            A.setEntry(2, 1, 0.0);
            A.setEntry(2, 2, 1.0);
            r.setEntry(2, 0.0);
        } else {
            // General case: n >= 4 points - use not-a-knot boundary conditions
            // Set up the system A*mu = r for the not-a-knot condition
            // First row: not-a-knot condition at x[1]
            A.setEntry(0, 0, h[1]);
            A.setEntry(0, 1, -(h[0] + h[1]));
            A.setEntry(0, 2, h[0]);
            
            // Last row: not-a-knot condition at x[n-2]
            A.setEntry(n - 1, n - 3, h[n - 2]);
            A.setEntry(n - 1, n - 2, -(h[n - 2] + h[n - 3]));
            A.setEntry(n - 1, n - 1, h[n - 3]);

            for (int i = 1; i < n - 1; i++) {
                A.setEntry(i, i - 1, h[i - 1]);
                A.setEntry(i, i, 2 * (h[i - 1] + h[i]));
                A.setEntry(i, i + 1, h[i]);
                r.setEntry(i, 6 * (delta[i] - delta[i - 1]));
            }
        }

        // Solve for mu
        DecompositionSolver solver = new LUDecomposition(A).getSolver();
        RealVector mu = solver.solve(r);

        // Store coefficients for each spline segment
        // logI_i(t) = splineCoeffs[i][3] + splineCoeffs[i][2]*(t-x[i]) + splineCoeffs[i][1]*(t-x[i])^2 + splineCoeffs[i][0]*(t-x[i])^3
        splineCoeffs = new double[n - 1][4];
        for (int i = 0; i < n - 1; i++) {
            splineCoeffs[i][0] = (mu.getEntry(i + 1) - mu.getEntry(i)) / (6 * h[i]);
            splineCoeffs[i][1] = mu.getEntry(i) / 2;
            splineCoeffs[i][2] = delta[i] - h[i] * (2 * mu.getEntry(i) + mu.getEntry(i + 1)) / 6;
            splineCoeffs[i][3] = infected.getArrayValue(i);
        }
    }
}