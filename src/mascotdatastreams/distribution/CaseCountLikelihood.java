package mascotdatastreams.distribution;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.parameter.RealParameter;
import mascot.parameterdynamics.Skygrowth;

import java.util.List;
import java.util.Random;

public class CaseCountLikelihood extends Distribution {
    // Single-deme prevalence input provided via Mascot's Skygrowth.
    // Note: Skygrowth#getNeTime(t) is interpreted here as I(t), i.e., prevalence at time t in forward time.
    public final Input<Skygrowth> prevalenceSingleInput = new Input<>(
            "prevalence", "Mascot Skygrowth providing prevalence values I(t) via getNeTime(t)", Validate.REQUIRED);

    // Observations passed directly as parameters
    public final Input<RealParameter> caseCountsInput = new Input<>(
            "caseCounts", "Observed case counts (dimension = number of observations)", Validate.REQUIRED);
    public final Input<RealParameter> caseTimesInput = new Input<>(
            "caseTimes", "Observation times corresponding 1:1 to caseCounts", Validate.REQUIRED);
    final public Input<ParametricDistribution> distInput = new Input<>("distribution", "Distribution used to calculate likelihood. Currently only GammaPoisson is supported.", Validate.REQUIRED);
    // Optional placeholder if needed elsewhere; not used directly in this likelihood
    public final Input<RealParameter> uninfectiousRateInput = new Input<>(
            "uninfectiousRate", "Fixed uninfectious rate (per time unit), optional.", Validate.OPTIONAL);
    // Optional scaling of the mean (akin to a sampling/surveillance rate in a given deme); defaults to 1.0 (no scaling).
    public final Input<RealParameter> scalingInput = new Input<>(
            "scaling", "Scaling factor applied to the prevalence-derived mean; must be > 0.", Validate.OPTIONAL);
    
    protected RealParameter caseCounts;
    protected RealParameter caseTimes;
    protected ParametricDistribution dist;
    protected Skygrowth prevalence;
    protected boolean validated = false;
    
    @Override
    public void initAndValidate() {
        caseCounts = caseCountsInput.get();
        caseTimes = caseTimesInput.get();
        dist = distInput.get();
        prevalence = prevalenceSingleInput.get();

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
        for (int i = 0; i < n; i++) {
            double t = caseTimes.getArrayValue(i);
            double caseCount = caseCounts.getArrayValue(i);

            // Get prevalence at this time via Skygrowth and take log for clarity
            // getNeTime provides I(t) under our prevalence semantics
            double meanI = prevalence.getNeTime(t);

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
            if (!(dist instanceof CountDistributionWithMean)) {
                throw new IllegalArgumentException(
                        "CaseCountLikelihood requires distributions implementing CountDistributionWithMean. Got: "
                                + dist.getClass().getName());
            }
            double intervalLogP = ((CountDistributionWithMean) dist).logPForMean(caseCount, scaledMean);
            if (Double.isNaN(intervalLogP)) {
                System.err.println("Warning: NaN likelihood for observation index " + i + ": caseCount=" + caseCount + ", scaled mean=" + scaledMean);
                return Double.NEGATIVE_INFINITY;
            }

            logP += intervalLogP;
        }
        return logP;
    }

    @Override
    public boolean requiresRecalculation() {
        boolean dirty = false;
        // Upstream process (Skygrowth) driving prevalence
        if (prevalence != null && prevalence.isDirty()) dirty = true;
        // Observations
        if (caseCounts != null) {
            int n = caseCounts.getDimension();
            for (int i = 0; i < n; i++) {
                if (caseCounts.isDirty(i)) { dirty = true; break; }
            }
        }
        if (caseTimes != null) {
            int n = caseTimes.getDimension();
            for (int i = 0; i < n; i++) {
                if (caseTimes.isDirty(i)) { dirty = true; break; }
            }
        }
        // Optional parameters
        RealParameter scaleParam = scalingInput.get();
        if (scaleParam != null && scaleParam.isDirty(0)) dirty = true;
        RealParameter uninf = uninfectiousRateInput.get();
        if (uninf != null && uninf.isDirty(0)) dirty = true;
        // Distribution's dispersion parameter (if GammaPoisson with RealParameter dispersion)
        if (dist instanceof GammaPoisson) {
            GammaPoisson gp = (GammaPoisson) dist;
            if (gp.dispersionInput.get() instanceof RealParameter) {
                RealParameter disp = (RealParameter) gp.dispersionInput.get();
                if (disp != null && disp.isDirty(0)) dirty = true;
            }
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
