package mascotdatastreams.dynamics;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.RealParameter;
import mascot.dynamics.RateShifts;
import org.apache.commons.math3.linear.*;


/**
 * @author Nicola F. Mueller
 */
@Description("Populaiton function with values at certain time points that are interpolated in between. Parameter has to be in log space")
public class NotAKnotSpline extends CalculationNode {

    final public Input<RealParameter> InfectedInput = new Input<>("logInfected",
            "Nes over time in log space", Input.Validate.REQUIRED);
    final public Input<RateShifts> rateShiftsInput = new Input<>("rateShifts",
            "When to switch between elements of Ne", Input.Validate.REQUIRED);
    final public Input<RateShifts> gridRateShiftsInput = new Input<>("gridRateShifts",
            "Rate shifts to use as grid points for spline evaluation (optional, defaults to rateShifts)", Input.Validate.REQUIRED);
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
        infected.setDimension(rateShifts.getDimension());
        uninfectiousRate = uninfectiousRateInput.get();
        recalculateRates();
    }

    // computes the Ne's at the break points from the growth rates and the transmission rates
    private void recalculateRates() {
        notAKnotCubicSpline();
        // use grid rate shifts as grid points
        int n = gridRateShifts.getDimension() + 1; // +1 for time 0
        time = new double[n];
        I = new double[n];
        transmissionRate = new double[n];
        
        isValid = true;
        
        // First point at time 0
        time[0] = 0.0;
        I[0] = Math.exp(infected.getArrayValue(0));
        transmissionRate[0] = uninfectiousRate.getValue() - 
                I[0] * splineCoeffs[0][2]; // derivative at t=0 is splineCoeffs[0][2]
        if (transmissionRate[0] < 0) {
            isValid = false;
        }
        
        // Points at grid rate shifts
        for (int i = 0; i < gridRateShifts.getDimension(); i++) {
            
            time[i+1] = gridRateShifts.getValue(i);
            // Find which spline segment this grid point falls into
            int segmentIndex = findSplineSegment(gridRateShifts.getValue(i));
            
            // Bounds check for segmentIndex
            if (segmentIndex >= splineCoeffs.length) {
                segmentIndex = splineCoeffs.length - 1;
            }
                        
            double timeDiff = gridRateShifts.getValue(i) - rateShifts.getValue(segmentIndex);
            double timeDiff2 = timeDiff * timeDiff;
            double timeDiff3 = timeDiff2 * timeDiff;
            
            // Evaluate spline at this grid point
            I[i+1] = Math.exp(splineCoeffs[segmentIndex][0]*timeDiff3 + splineCoeffs[segmentIndex][1]*timeDiff2 + 
                             splineCoeffs[segmentIndex][2]*timeDiff + splineCoeffs[segmentIndex][3]);
            
            // Calculate derivative at this grid point
            double derivative = 3 * splineCoeffs[segmentIndex][0] * timeDiff2 + 2 * splineCoeffs[segmentIndex][1] * timeDiff + splineCoeffs[segmentIndex][2];
            transmissionRate[i+1] = uninfectiousRate.getValue() - I[i+1] * derivative;
            
            if (transmissionRate[i+1] < 0) {
                isValid = false;
            }
        }
        
        ratesKnows = true;
    }
    
    // Helper method to find which spline segment a given time falls into
    private int findSplineSegment(double t) {
        for (int i = 0; i < rateShifts.getDimension(); i++) {
            if (t < rateShifts.getValue(i)) {
                return i;
            }
        }
        // Use last segment if t is beyond all rate shifts
        // splineCoeffs has rateShifts.getDimension() - 1 elements (indices 0 to rateShifts.getDimension()-2)
        return rateShifts.getDimension() - 2;
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
            int mid = (left + right + 1) / 2;
            if (time[mid] <= t) {
                left = mid;
            } else {
                right = mid - 1;
            }
        }
        
        return Math.log(I[left]);
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
            int k = 0;
            double timeDiff = t - rateShifts.getValue(k);
            double timeDiff2 = timeDiff * timeDiff;
            return 3 * splineCoeffs[k][0] * timeDiff2 + 2 * splineCoeffs[k][1] * timeDiff + splineCoeffs[k][2];
        }
        
        // Binary search for the closest grid point
        int left = 0, right = time.length - 1;
        while (left < right) {
            int mid = (left + right + 1) / 2;
            if (time[mid] <= t) {
                left = mid;
            } else {
                right = mid - 1;
            }
        }
        
        // Calculate derivative using spline coefficients
        int k = left;
        if (k >= rateShifts.getDimension() - 1) {
            k = rateShifts.getDimension() - 2;
        }
        double timeDiff = t - rateShifts.getValue(k);
        double timeDiff2 = timeDiff * timeDiff;
        return 3 * splineCoeffs[k][0] * timeDiff2 + 2 * splineCoeffs[k][1] * timeDiff + splineCoeffs[k][2];
    }
    
    /**
     * Get the number of grid points in the spline.
     * 
     * @return number of grid points
     */
    public int getGridPointCount() {
        return gridRateShifts.getDimension() + 1; // +1 for time 0
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
    
    /**
     * Get the step size between grid points.
     * Since grid points are now at grid rate shifts, this returns the average step size.
     * 
     * @return average step size
     */
    public double getGridStep() {
        if (gridRateShifts.getDimension() > 0) {
            return getGridEnd() / gridRateShifts.getDimension();
        }
        return 0.0;
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

        // Calculate h values (difference between x values)
        double[] h = new double[n - 1];
        for (int i = 0; i < n - 1; i++) {
            h[i] = rateShifts.getValue(i + 1) - rateShifts.getValue(i);
        }

        // Calculate the difference in y values
        double[] delta = new double[n - 1];
        for (int i = 0; i < n - 1; i++) {
            delta[i] = (infected.getArrayValue(i + 1) - infected.getArrayValue(i)) / h[i];
        }

        // Create the tridiagonal system
        RealMatrix A = new Array2DRowRealMatrix(n, n);
        RealVector r = new ArrayRealVector(n);

        // Set up the system A*mu = r for the not-a-knot condition
        A.setEntry(0, 0, h[1]);
        A.setEntry(0, 1, -(h[0] + h[1]));
        A.setEntry(0, 2, h[0]);
        A.setEntry(n - 1, n - 3, h[n - 2]);
        A.setEntry(n - 1, n - 2, -(h[n - 2] + h[n - 3]));
        A.setEntry(n - 1, n - 1, h[n - 3]);

        for (int i = 1; i < n - 1; i++) {
            A.setEntry(i, i - 1, h[i - 1]);
            A.setEntry(i, i, 2 * (h[i - 1] + h[i]));
            A.setEntry(i, i + 1, h[i]);
            r.setEntry(i, 6 * (delta[i] - delta[i - 1]));
        }

        // Solve for mu
        DecompositionSolver solver = new LUDecomposition(A).getSolver();
        RealVector mu = solver.solve(r);

        // Store coefficients for each spline segment
        splineCoeffs = new double[n - 1][4];
        for (int i = 0; i < n - 1; i++) {
            splineCoeffs[i][0] = (mu.getEntry(i + 1) - mu.getEntry(i)) / (6 * h[i]);
            splineCoeffs[i][1] = mu.getEntry(i) / 2;
            splineCoeffs[i][2] = delta[i] - h[i] * (2 * mu.getEntry(i) + mu.getEntry(i + 1)) / 6;
            splineCoeffs[i][3] = infected.getArrayValue(i);
        }
    }
}