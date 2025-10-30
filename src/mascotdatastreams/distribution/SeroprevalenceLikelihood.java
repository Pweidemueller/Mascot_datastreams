package mascotdatastreams.distribution;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.parameter.RealParameter;
import mascotdatastreams.dynamics.Spline;

import java.util.List;
import java.util.Random;

/**
 * Likelihood for cumulative case count observations based on prevalence dynamics from a spline.
 * 
 * This class takes a Spline that provides log-prevalence values and computes
 * the cumulative cases by integrating the prevalence over time. The cumulative cases
 * are then used to calculate the likelihood of observed cumulative case counts.
 * 
 * The integral calculation follows the template pattern where:
 * - Prevalence I(t) is obtained from the spline
 * - Cumulative cases = integral of I(t) from start time to observation time
 * - The likelihood is computed using the cumulative cases as the mean
 */
public class SeroprevalenceLikelihood extends Distribution {
    
    // Spline input providing log-prevalence values
    public final Input<Spline> prevalenceSplineInput = new Input<>(
            "prevalenceSpline", "Spline providing log-prevalence values for interpolation", Validate.REQUIRED);

    // Observations passed directly as parameters
    public final Input<RealParameter> SeroPeopleTestedInput = new Input<>(
            "seroPeopleTested", "Observed number of people tested for seropositivity (dimension = number of observations)", Validate.REQUIRED);
    public final Input<RealParameter> SeroPeopleSeropositiveInput = new Input<>(
            "seroPeopleSeropositive", "Observed number of people seropositive (dimension = number of observations)", Validate.REQUIRED);
    public final Input<RealParameter> SeroTimesInput = new Input<>(
            "seroTimes", "Observation times corresponding 1:1 to seroPeopleSeropositive", Validate.REQUIRED);
    
    final public Input<ParametricDistribution> distInput = new Input<>(
            "distribution", "Distribution used to calculate likelihood. Expected Binomial.", Validate.REQUIRED);
    
    // Optional log-space scaling factor
    public final Input<RealParameter> logScalingInput = new Input<>(
            "scaling", "Scaling factor applied to the prevalence-derived mean in LOG space.", Validate.OPTIONAL);
    
    // Start time for cumulative calculation (default: earliest time in spline)
    public final Input<RealParameter> startTimeInput = new Input<>(
            "startTime", "Start time for cumulative case calculation (default: earliest spline time)", Validate.OPTIONAL);
    
    protected RealParameter SeroPeopleTested;
    protected RealParameter SeroPeopleSeropositive;
    protected RealParameter SeroTimes;
    protected RealParameter logScaling;
    protected ParametricDistribution dist;
    protected Spline prevalenceSpline;

    protected boolean validated = false;
    
    @Override
    public void initAndValidate() {
        SeroPeopleTested = SeroPeopleTestedInput.get();
        SeroPeopleSeropositive = SeroPeopleSeropositiveInput.get();
        SeroTimes = SeroTimesInput.get();
        dist = distInput.get();
        prevalenceSpline = prevalenceSplineInput.get();

        if (!(dist instanceof Binomial)) {
            throw new IllegalArgumentException("SeroprevalenceLikelihood: 'distribution' must be Binomial");
        }

        // Validate optional scaling parameter if present
        logScaling = logScalingInput.get();
        if (logScaling != null) {
            double s = Math.exp(logScaling.getArrayValue());
            if (!(s > 0.0)) {
                throw new IllegalArgumentException("SeroprevalenceLikelihood: 'logScaling' must be > 0, got " + s);
            }
        }

        // Validate observation vectors
        if (SeroPeopleTested == null || SeroPeopleSeropositive == null || SeroTimes == null) {
            throw new IllegalArgumentException("SeroprevalenceLikelihood: 'SeroPeopleTested', 'SeroPeopleSeropositive' and 'SeroTimes' must be provided");
        }
        if (SeroPeopleTested.getDimension() != SeroPeopleSeropositive.getDimension() || SeroPeopleTested.getDimension() != SeroTimes.getDimension()) {
            throw new IllegalArgumentException("SeroprevalenceLikelihood: 'SeroPeopleTested', 'SeroPeopleSeropositive' and 'SeroTimes' must have the same dimension");
        }
        if (SeroPeopleTested.getDimension() == 0) {
            System.err.println("Warning: SeroprevalenceLikelihood received zero observations. Likelihood will be neutral (0).");
        }

        calculateLogP();
        validated = true;
    }
    
    /**
     * Computes the integral of prevalence from startTime to endTime.
     * This represents the cumulative cases over the time period.
     * 
     * @param startTime start time for integration
     * @param endTime end time for integration
     * @return cumulative cases (integral of prevalence)
     */
    private double getCumulativeCases(double startTime, double endTime) {
        if (startTime >= endTime) {
            return 0.0;
        }
        
        // Grid info
        int n = prevalenceSpline.getGridPointCount();
        if (n <= 1) {
            return 0.0;
        }

        double gridMin = prevalenceSpline.getGridStart();
        double gridMax = prevalenceSpline.getGridEnd();

        // Clamp to grid range
        double from = Math.max(startTime, gridMin);
        double to = Math.min(endTime, gridMax);
        if (from >= to) {
            return 0.0;
        }

        // Find boundary indices
        int leftIdx = prevalenceSpline.getLeftGridIndex(from);
        int rightIdx = prevalenceSpline.getRightGridIndex(to);
        // Integrate using trapezoids over gridpoint segments
        double sum = 0.0;

        for (int i = leftIdx; i < rightIdx; i++) {
            double t1 = prevalenceSpline.getGridPointTime(i);
            double t2 = prevalenceSpline.getGridPointTime(i + 1);
            double I1 = prevalenceSpline.getValueAtGridPoint(t1);
            double I2 = prevalenceSpline.getValueAtGridPoint(t2);

            sum += (t2 - t1) * (I1 + I2) * 0.5;
        }

        return sum;
    }
    
    @Override
    public double calculateLogP() {
        logP = 0.0;
        // Neutral likelihood if there are no observations
        if (SeroPeopleTested == null || SeroPeopleTested.getDimension() == 0) {
            return logP;
        }

        double startTime = prevalenceSpline.getGridStart(); // Default to earliest grid time     
        RealParameter startTimeParam = startTimeInput.get();
        if (startTimeParam != null) {
            startTime = startTimeParam.getArrayValue();
        }

        int numObservations = SeroPeopleTested.getDimension();
        for (int i = 0; i < numObservations; i++) {
            double t = SeroTimes.getArrayValue(i);
            int n = (int) Math.round(SeroPeopleTested.getArrayValue(i));
            int x = (int) Math.round(SeroPeopleSeropositive.getArrayValue(i));

            // Determine start time for cumulative calculation


            // Get cumulative cases from start time to observation time
            double cumulativeCases = getCumulativeCases(startTime, t);

            // Apply optional scaling factor (default 1.0)
            double scaling = 1.0;
            if (logScaling != null) {
                scaling = Math.exp(logScaling.getArrayValue());
                if (!(scaling > 0.0)) {
                    System.err.println("Warning: Invalid 'scaling' <= 0 encountered: " + scaling);
                    return Double.NEGATIVE_INFINITY;
                }
            }
            // Convert cumulative cases (expected count in whole population) to probability via scaling
            double p = cumulativeCases * scaling;
            // Clamp to (0,1) for numerical stability and validity
            final double PROB_EPS = 1e-16;
            p = Math.min(1.0 - PROB_EPS, Math.max(PROB_EPS, p));

            // Validate observation bounds
            if (n < 0 || x < 0 || x > n) {
                System.err.println("Warning: Invalid observation at index " + i + ": n=" + n + ", x=" + x);
                return Double.NEGATIVE_INFINITY;
            }

            // Calculate binomial log-likelihood using the provided distribution instance
            double intervalLogP = ((Binomial) dist).logPMFForParams(x, n, p);
            if (Double.isNaN(intervalLogP)) {
                System.err.println("Warning: NaN likelihood for observation index " + i + ": n=" + n + ", x=" + x + ", p=" + p);
                return Double.NEGATIVE_INFINITY;
            }

            logP += intervalLogP;
        }
        return logP;
    }

    @Override
    public boolean requiresRecalculation() {
        return true;
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
        // Not implemented for likelihood classes
    }
}
