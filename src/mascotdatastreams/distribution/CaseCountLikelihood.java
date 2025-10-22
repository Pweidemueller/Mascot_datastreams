package mascotdatastreams.distribution;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import mascotdatastreams.dynamics.NotAKnotSpline;
import org.apache.commons.math.special.Gamma;

import java.util.List;
import java.util.Random;

/**
 * Likelihood for case count observations based on prevalence dynamics from a spline.
 * 
 * This class takes a NotAKnotSpline that provides log-prevalence values and converts
 * them to prevalence I(t) for calculating the likelihood of observed case counts.
 * The spline interpolation allows for smooth prevalence trajectories between
 * rate shift points.
 */
public class CaseCountLikelihood extends Distribution {
    // Single-deme prevalence input provided via NotAKnotSpline.
    // Note: The spline provides log-prevalence values that are converted to prevalence I(t).
    public final Input<NotAKnotSpline> prevalenceSplineInput = new Input<>(
            "prevalenceSpline", "NotAKnotSpline providing log-prevalence values for interpolation", Validate.REQUIRED);

    // Observations passed directly as parameters
    public final Input<RealParameter> caseCountsInput = new Input<>(
            "caseCounts", "Observed case counts (dimension = number of observations)", Validate.REQUIRED);
    public final Input<RealParameter> caseTimesInput = new Input<>(
            "caseTimes", "Observation times corresponding 1:1 to caseCounts", Validate.REQUIRED);
    // Dispersion parameter for Gamma-Poisson (Negative Binomial) distribution
    public final Input<RealParameter> dispersionInput = new Input<>(
            "dispersion", "Dispersion parameter (alpha) for Gamma-Poisson distribution", Validate.REQUIRED);
    // Optional placeholder if needed elsewhere; not used directly in this likelihood
    public final Input<RealParameter> uninfectiousRateInput = new Input<>(
            "uninfectiousRate", "Fixed uninfectious rate (per time unit), optional.", Validate.OPTIONAL);
    // Optional scaling of the mean (akin to a sampling/surveillance rate in a given deme); defaults to 1.0 (no scaling).
    public final Input<RealParameter> scalingInput = new Input<>(
            "scaling", "Scaling factor applied to the prevalence-derived mean; must be > 0.", Validate.OPTIONAL);
    
    protected RealParameter caseCounts;
    protected RealParameter caseTimes;
    protected RealParameter dispersion;
    protected NotAKnotSpline prevalenceSpline;
    protected boolean validated = false;
    
    @Override
    public void initAndValidate() {
        caseCounts = caseCountsInput.get();
        caseTimes = caseTimesInput.get();
        dispersion = dispersionInput.get();
        prevalenceSpline = prevalenceSplineInput.get();

        // Validate dispersion parameter
        if (dispersion == null) {
            throw new IllegalArgumentException("CaseCountLikelihood: 'dispersion' parameter is required");
        }
        double alpha = dispersion.getArrayValue();
        if (!(alpha > 0.0)) {
            throw new IllegalArgumentException("CaseCountLikelihood: 'dispersion' must be > 0, got " + alpha);
        }

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
    
    @Override
    public double calculateLogP() {
        logP = 0.0;
        // Neutral likelihood if there are no observations
        if (caseCounts == null || caseCounts.getDimension() == 0) {
            return logP;
        }

        int n = caseCounts.getDimension();
//        System.out.println("Calculating CaseCountLikelihood for " + n + " observations.");
        for (int i = 0; i < n; i++) {
            double t = caseTimes.getArrayValue(i);
            double caseCount = caseCounts.getArrayValue(i);

            // Get prevalence at this time via spline interpolation
            // The spline provides log-prevalence values, so we need to exponentiate
            double logI = prevalenceSpline.getValueAtGridPoint(t);
            double meanI = Math.exp(logI);
//            System.out.println(meanI + " " + caseCount);

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


            // Calculate Gamma-Poisson (Negative Binomial) log likelihood directly
            double alpha = dispersion.getArrayValue();
            double intervalLogP = calculateGammaPoissonLogP(caseCount, scaledMean, alpha);
            if (Double.isNaN(intervalLogP)) {
                System.err.println("Warning: NaN likelihood for observation index " + i + ": caseCount=" + caseCount + ", scaled mean=" + scaledMean + ", alpha=" + alpha);
                return Double.NEGATIVE_INFINITY;
            }

            logP += intervalLogP;
        }
        return logP;
    }

    /**
     * Calculate the log probability of a Gamma-Poisson (Negative Binomial) distribution.
     * 
     * @param observation observed count
     * @param mean mean of the distribution
     * @param alpha dispersion parameter
     * @return log probability
     */
    private double calculateGammaPoissonLogP(double observation, double mean, double alpha) {
        
        int x = (int) Math.round(observation);
        double r = 1.0 / alpha;
        double p = r / (r + mean);
        
        // Clamp p to avoid numerical issues
        p = Math.min(1 - 1e-16, Math.max(1e-16, p));
        
        // Log PMF = ln Γ(r + x) - ln Γ(r) - ln Γ(x+1) + r ln p + x ln(1-p)
        double logGammaRK = Gamma.logGamma(r + x);
        double logGammaR = Gamma.logGamma(r);
        double logGammaK1 = Gamma.logGamma(x + 1.0);
        double logP = Math.log(p);
        double log1mP = Math.log(1.0 - p);
        
        return logGammaRK - logGammaR - logGammaK1 + r * logP + x * log1mP;
    }

    @Override
    public boolean requiresRecalculation() {
        boolean dirty = false;
        
        // Check observations
        if (caseCounts != null) {
            int n = caseCounts.getDimension();
            for (int i = 0; i < n; i++) {
                if (caseCounts.isDirty(i)) { 
                    dirty = true; 
                    break; 
                }
            }
        }
        if (caseTimes != null) {
            int n = caseTimes.getDimension();
            for (int i = 0; i < n; i++) {
                if (caseTimes.isDirty(i)) { 
                    dirty = true; 
                    break; 
                }
            }
        }
        
        // Check dispersion parameter
        if (dispersion != null && dispersion.isDirty(0)) {
            dirty = true;
        }
        
        // Check optional parameters
        RealParameter scaleParam = scalingInput.get();
        if (scaleParam != null && scaleParam.isDirty(0)) {
            dirty = true;
        }
        
        return dirty || super.requiresRecalculation();
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
