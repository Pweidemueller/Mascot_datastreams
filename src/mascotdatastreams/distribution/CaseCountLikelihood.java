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
 * Likelihood for case count observations based on prevalence dynamics from a spline.
 * 
 * This class takes a Spline that provides log-prevalence values and converts
 * them to prevalence I(t) for calculating the likelihood of observed case counts.
 * The spline interpolation allows for smooth prevalence trajectories between
 * rate shift points.
 */
public class CaseCountLikelihood extends Distribution {
    // Single-deme prevalence input provided via Spline.
    // Note: The spline provides log-prevalence values that are converted to prevalence I(t).
    public final Input<Spline> prevalenceSplineInput = new Input<>(
            "prevalenceSpline", "Spline providing log-prevalence values for interpolation", Validate.REQUIRED);

    // Observations passed directly as parameters
    public final Input<RealParameter> caseCountsInput = new Input<>(
            "caseCounts", "Observed case counts (dimension = number of observations)", Validate.REQUIRED);
    public final Input<RealParameter> caseTimesInput = new Input<>(
            "caseTimes", "Observation times corresponding 1:1 to caseCounts", Validate.REQUIRED);
    final public Input<ParametricDistribution> distInput = new Input<>("distribution", "Distribution used to calculate likelihood. Currently only GammaPoisson is supported.", Validate.REQUIRED);
    // Optional scaling of the mean (akin to a sampling/surveillance rate in a given deme); defaults to 1.0 (no scaling).
    public final Input<RealParameter> scalingInput = new Input<>(
            "scaling", "Scaling factor applied to the prevalence-derived mean; must be > 0.", Validate.OPTIONAL);
    
    protected RealParameter caseCounts;
    protected RealParameter caseTimes;
    protected ParametricDistribution dist;
    protected Spline prevalenceSpline;
    protected boolean validated = false;
    
    @Override
    public void initAndValidate() {
        caseCounts = caseCountsInput.get();
        caseTimes = caseTimesInput.get();
        dist = distInput.get();
        prevalenceSpline = prevalenceSplineInput.get();

        // Validate optional scaling parameter if present
        RealParameter scaleParam = scalingInput.get();
        if (scaleParam != null) {
            double s = scaleParam.getArrayValue();
            if (!(s > 0.0)) {
                throw new IllegalArgumentException("CaseCountLikelihood: 'scaling' must be > 0, got " + s);
            }
        }

        // Validate observation vectors
        if (caseCounts == null || caseTimes == null) {
            throw new IllegalArgumentException("CaseCountLikelihood: 'caseCounts' and 'caseTimes' must be provided");
        }
        if (caseCounts.getDimension() != caseTimes.getDimension()) {
            throw new IllegalArgumentException("CaseCountLikelihood: 'caseCounts' and 'caseTimes' must have the same dimension");
        }
        if (caseCounts.getDimension() == 0) {
            System.err.println("Warning: CaseCountLikelihood received zero observations. Likelihood will be neutral (0).");
        }

        calculateLogP();
        validated = true;
    }
    
    int ii=0;
    
    @Override
    public double calculateLogP() {
    	ii++;
        logP = 0.0;
        // Neutral likelihood if there are no observations
        if (caseCounts == null || caseCounts.getDimension() == 0) {
            return logP;
        }

        int n = caseCounts.getDimension();
        for (int i = 0; i < n; i++) {
            double t = caseTimes.getArrayValue(i);
            double caseCount = caseCounts.getArrayValue(i);

            // Get prevalence at this time via spline interpolation
            // The spline provides log-prevalence values, so we need to exponentiate
            double logI = prevalenceSpline.getValueAtGridPoint(t);
            double meanI = Math.exp(logI);

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
            double scaledMean = meanI * scaling;

            // Validate parameters
            if (scaledMean <= 0.0 || caseCount < 0.0) {
                System.err.println("Warning: Invalid parameters - scaled mean: " + scaledMean + ", caseCount: " + caseCount);
                return Double.NEGATIVE_INFINITY;
            }

            // Calculate log likelihood using a stateless interface (no input mutation)
            if (!(dist instanceof DistributionWithMean)) {
                throw new IllegalArgumentException(
                        "CaseCountLikelihood requires distributions implementing DistributionWithMean. Got: "
                                + dist.getClass().getName());
            }
            double intervalLogP = ((DistributionWithMean) dist).logPForMean(caseCount, scaledMean);
            if (Double.isNaN(intervalLogP)) {
                System.err.println("Warning: NaN likelihood for observation index " + i + ": caseCount=" + caseCount + ", scaled mean=" + scaledMean);
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
