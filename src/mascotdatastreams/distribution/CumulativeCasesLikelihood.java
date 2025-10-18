package mascotdatastreams.distribution;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.parameter.RealParameter;
import mascotdatastreams.dynamics.NotAKnotSpline;

import java.util.List;
import java.util.Random;

/**
 * Likelihood for cumulative case count observations based on prevalence dynamics from a spline.
 * 
 * This class takes a NotAKnotSpline that provides log-prevalence values and computes
 * the cumulative cases by integrating the prevalence over time. The cumulative cases
 * are then used to calculate the likelihood of observed cumulative case counts.
 * 
 * The integral calculation follows the template pattern where:
 * - Prevalence I(t) is obtained from the spline
 * - Cumulative cases = integral of I(t) from start time to observation time
 * - The likelihood is computed using the cumulative cases as the mean
 */
public class CumulativeCasesLikelihood extends Distribution {
    
    // Spline input providing log-prevalence values
    public final Input<NotAKnotSpline> prevalenceSplineInput = new Input<>(
            "prevalenceSpline", "NotAKnotSpline providing log-prevalence values for interpolation", Validate.REQUIRED);

    // Observations passed directly as parameters
    public final Input<RealParameter> cumulativeCaseCountsInput = new Input<>(
            "cumulativeCaseCounts", "Observed cumulative case counts (dimension = number of observations)", Validate.REQUIRED);
    public final Input<RealParameter> caseTimesInput = new Input<>(
            "caseTimes", "Observation times corresponding 1:1 to cumulativeCaseCounts", Validate.REQUIRED);
    
    // Distribution for likelihood calculation
    final public Input<ParametricDistribution> distInput = new Input<>("distribution", 
            "Distribution used to calculate likelihood. Currently only GammaPoisson is supported.", Validate.REQUIRED);
    
    // Optional parameters
    public final Input<RealParameter> uninfectiousRateInput = new Input<>(
            "uninfectiousRate", "Fixed uninfectious rate (per time unit), optional.", Validate.OPTIONAL);
    public final Input<RealParameter> scalingInput = new Input<>(
            "scaling", "Scaling factor applied to the prevalence-derived mean; must be > 0.", Validate.OPTIONAL);
    
    // Start time for cumulative calculation (default: earliest time in spline)
    public final Input<RealParameter> startTimeInput = new Input<>(
            "startTime", "Start time for cumulative case calculation (default: earliest spline time)", Validate.OPTIONAL);
    
    protected RealParameter cumulativeCaseCounts;
    protected RealParameter caseTimes;
    protected ParametricDistribution dist;
    protected NotAKnotSpline prevalenceSpline;
    protected boolean validated = false;
    
    @Override
    public void initAndValidate() {
        cumulativeCaseCounts = cumulativeCaseCountsInput.get();
        caseTimes = caseTimesInput.get();
        dist = distInput.get();
        prevalenceSpline = prevalenceSplineInput.get();

        // Validate optional scaling parameter if present
        RealParameter scaleParam = scalingInput.get();
        if (scaleParam != null) {
            double s = scaleParam.getArrayValue();
            if (!(s > 0.0)) {
                throw new IllegalArgumentException("CumulativeCasesLikelihood: 'scaling' must be > 0, got " + s);
            }
        }

        // Validate observation vectors
        if (cumulativeCaseCounts == null || caseTimes == null) {
            throw new IllegalArgumentException("CumulativeCasesLikelihood: 'cumulativeCaseCounts' and 'caseTimes' must be provided");
        }
        if (cumulativeCaseCounts.getDimension() != caseTimes.getDimension()) {
            throw new IllegalArgumentException("CumulativeCasesLikelihood: 'cumulativeCaseCounts' and 'caseTimes' must have the same dimension");
        }
        if (cumulativeCaseCounts.getDimension() == 0) {
            System.err.println("Warning: CumulativeCasesLikelihood received zero observations. Likelihood will be neutral (0).");
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
        
        // Get spline grid information
        int gridPoints = prevalenceSpline.getGridPointCount();
        double gridStart = prevalenceSpline.getGridStart();
        double gridEnd = prevalenceSpline.getGridEnd();
        double gridStep = prevalenceSpline.getGridStep();
        
        // Clamp times to spline range
        double from = Math.max(startTime, gridStart);
        double to = Math.min(endTime, gridEnd);
        
        if (from >= to) {
            return 0.0;
        }
        
        // Compute the integral using trapezoidal rule
        double cumulativeCases = 0.0;
        
        // Find the starting grid point
        int startIndex = (int) Math.max(0, Math.floor((from - gridStart) / gridStep));
        int endIndex = (int) Math.min(gridPoints - 1, Math.ceil((to - gridStart) / gridStep));
        
        for (int i = startIndex; i < endIndex; i++) {
            double t1 = gridStart + i * gridStep;
            double t2 = gridStart + (i + 1) * gridStep;
            
            // Clamp to integration bounds
            t1 = Math.max(t1, from);
            t2 = Math.min(t2, to);
            
            if (t1 >= t2) continue;
            
            // Get prevalence values at grid points
            double logI1 = prevalenceSpline.getValueAtGridPoint(t1);
            double logI2 = prevalenceSpline.getValueAtGridPoint(t2);
            double I1 = Math.exp(logI1);
            double I2 = Math.exp(logI2);
            
            // Trapezoidal rule: integral ≈ (t2-t1) * (I1 + I2) / 2
            cumulativeCases += (t2 - t1) * (I1 + I2) / 2.0;
        }
        
        return cumulativeCases;
    }
    
    @Override
    public double calculateLogP() {
        logP = 0.0;
        // Neutral likelihood if there are no observations
        if (cumulativeCaseCounts == null || cumulativeCaseCounts.getDimension() == 0) {
            return logP;
        }

        int n = cumulativeCaseCounts.getDimension();
        for (int i = 0; i < n; i++) {
            double t = caseTimes.getArrayValue(i);
            double cumulativeCaseCount = cumulativeCaseCounts.getArrayValue(i);

            // Determine start time for cumulative calculation
            double startTime = prevalenceSpline.getGridStart(); // Default to spline start
            RealParameter startTimeParam = startTimeInput.get();
            if (startTimeParam != null) {
                startTime = startTimeParam.getArrayValue();
            }

            // Get cumulative cases from start time to observation time
            double cumulativeCases = getCumulativeCases(startTime, t);

            // Apply optional scaling factor (default 1.0)
            double scaling = 1.0;
            RealParameter scaleParam = scalingInput.get();
            if (scaleParam != null) {
                scaling = scaleParam.getArrayValue();
                if (!(scaling > 0.0)) {
                    System.err.println("Warning: Invalid 'scaling' <= 0 encountered: " + scaling);
                    return Double.NEGATIVE_INFINITY;
                }
            }
            double scaledMean = cumulativeCases * scaling;

            // Validate parameters
            if (scaledMean <= 0.0 || cumulativeCaseCount < 0.0) {
                System.err.println("Warning: Invalid parameters - scaled mean: " + scaledMean + ", cumulativeCaseCount: " + cumulativeCaseCount);
                return Double.NEGATIVE_INFINITY;
            }

            // Calculate log likelihood using a stateless interface (no input mutation)
            if (!(dist instanceof CountDistributionWithMean)) {
                throw new IllegalArgumentException(
                        "CumulativeCasesLikelihood requires distributions implementing CountDistributionWithMean. Got: "
                                + dist.getClass().getName());
            }
            double intervalLogP = ((CountDistributionWithMean) dist).logPForMean(cumulativeCaseCount, scaledMean);
            if (Double.isNaN(intervalLogP)) {
                System.err.println("Warning: NaN likelihood for observation index " + i + ": cumulativeCaseCount=" + cumulativeCaseCount + ", scaled mean=" + scaledMean);
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
