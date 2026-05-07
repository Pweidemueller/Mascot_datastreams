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
 * Likelihood for wastewater concentration observations based on prevalence dynamics from a spline.
 * 
 * This class takes a Spline that provides log-prevalence values and converts
 * them to prevalence I(t) for calculating the likelihood of observed wastewater concentrations.
 * The spline interpolation allows for smooth prevalence trajectories between
 * rate shift points.
 * 
 * The generative model assumes:
 *   log(concentration) ~ Normal(μ, σ²),  where μ = log(α · I(t) / N)
 *
 * Equivalently, concentration ~ LogNormal(μ, σ²) with median = α · I(t) / N.
 * α is an optional scaling factor, I(t) is the prevalence, and N is the population size.
 */
public class WastewaterLikelihood extends Distribution {
    // Single-deme prevalence input provided via Spline.
    // Note: The spline provides log-prevalence values that are converted to prevalence I(t).
    public final Input<Spline> prevalenceSplineInput = new Input<>(
            "prevalenceSpline", "Spline providing log-prevalence values for interpolation", Validate.REQUIRED);

    // Observations passed directly as parameters
    public final Input<RealParameter> concentrationsInput = new Input<>(
            "concentrations", "Observed wastewater concentrations (dimension = number of observations)", Validate.REQUIRED);
    public final Input<RealParameter> concentrationTimesInput = new Input<>(
            "concentrationTimes", "Observation times corresponding 1:1 to concentrations", Validate.REQUIRED);
    public final Input<RealParameter> populationSizeInput = new Input<>(
            "populationSize", "Population size", Validate.REQUIRED);
    final public Input<ParametricDistribution> distInput = new Input<>(
            "distribution", "Distribution used to calculate likelihood. Currently only LogNormal is supported.", Validate.REQUIRED);
    // Optional scaling of the mean (akin to a sampling/surveillance rate in a given deme); defaults to 1.0 (no scaling).
    public final Input<RealParameter> scalingInput = new Input<>(
            "scaling", "Scaling factor applied to the prevalence-derived mean; must be > 0.", Validate.OPTIONAL);
    
    protected RealParameter concentrations;
    protected RealParameter concentrationTimes;
    protected RealParameter populationSize;
    protected ParametricDistribution dist;
    protected Spline prevalenceSpline;
    protected boolean validated = false;
    
    @Override
    public void initAndValidate() {
        concentrations = concentrationsInput.get();
        concentrationTimes = concentrationTimesInput.get();
        populationSize = populationSizeInput.get();
        dist = distInput.get();
        prevalenceSpline = prevalenceSplineInput.get();

        // Validate that distribution is LogNormal
        if (!(dist instanceof LogNormal)) {
            throw new IllegalArgumentException(
                    "WastewaterLikelihood requires LogNormal distribution. Got: " + dist.getClass().getName());
        }

        // Validate optional scaling parameter if present
        RealParameter scaleParam = scalingInput.get();
        if (scaleParam != null) {
            double s = scaleParam.getArrayValue();
            if (!(s > 0.0)) {
                throw new IllegalArgumentException("WastewaterLikelihood: 'scaling' must be > 0, got " + s);
            }
        }

        // Validate observation vectors
        if (concentrations == null || concentrationTimes == null) {
            throw new IllegalArgumentException("WastewaterLikelihood: 'concentrations' and 'concentrationTimes' must be provided");
        }
        if (concentrations.getDimension() != concentrationTimes.getDimension()) {
            throw new IllegalArgumentException("WastewaterLikelihood: 'concentrations' and 'concentrationTimes' must have the same dimension");
        }
        if (concentrations.getDimension() == 0) {
            System.err.println("Warning: WastewaterLikelihood received zero observations. Likelihood will be neutral (0).");
        }

        calculateLogP();
        validated = true;
    }
    
    @Override
    public double calculateLogP() {
        logP = 0.0;
        // Neutral likelihood if there are no observations
        if (concentrations == null || concentrations.getDimension() == 0) {
            return logP;
        }

        int n = concentrations.getDimension();
        for (int i = 0; i < n; i++) {
            double t = concentrationTimes.getArrayValue(i);
            double concentration = concentrations.getArrayValue(i);

            // Get prevalence at this time via spline interpolation
            // The spline provides log-prevalence values, so we need to exponentiate
            double logI = prevalenceSpline.getLogPrevalence(t);
            double I = Math.exp(logI);
            double I_N = I / populationSize.getArrayValue();

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
            double scaledMedian = I_N * scaling;

            if (concentration <= 0.0) {
                System.err.println("Warning: Invalid concentration: " + concentration);
                return Double.NEGATIVE_INFINITY;
            }

            // Convert real-space median to log-space mean μ = log(median).
            // logPForMean on LogNormal expects μ directly (mean of the underlying normal in log space).
            double mu = Math.log(scaledMedian);

            if (!(dist instanceof DistributionWithMean)) {
                throw new IllegalArgumentException(
                        "WastewaterLikelihood requires distributions implementing DistributionWithMean. Got: "
                                + dist.getClass().getName());
            }
            double intervalLogP = ((DistributionWithMean) dist).logPForMean(concentration, mu);
            if (Double.isNaN(intervalLogP)) {
                System.err.println("Warning: NaN likelihood for observation index " + i + ": concentration=" + concentration + ", scaled mean=" + scaledMedian);
                return Double.NEGATIVE_INFINITY;
            }

            logP += intervalLogP;
        }
        return logP;
    }
    
    // TODO: think of a more granular decision on when to recalculate
    @Override
    public boolean requiresRecalculation() {
        return true;
    }

    @Override
    public List<String> getArguments() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public List<String> getConditions() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public void sample(State state, Random random) {
        // TODO Auto-generated method stub
        
    }
    
}
