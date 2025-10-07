package mascotdatastreams.distribution;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.parameter.RealParameter;
import mascotdatastreams.dynamics.PrevalenceSkygrowth;

import java.util.List;
import java.util.Random;
import java.util.ArrayList;

@Description("Calculates likelihood of observing case counts for a single deme given its prevalence trajectory using a specified distribution")
public class CaseCountLikelihood extends Distribution {
    // Single-deme prevalence input (explicitly PrevalenceSkygrowth)
    public final Input<PrevalenceSkygrowth> prevalenceSingleInput = new Input<>(
            "prevalence", "PrevalenceSkygrowth object for the selected deme", Validate.REQUIRED);

    final public Input<CaseCountData> caseCountInput = new Input<>("caseCounts", "Case count data", Validate.REQUIRED);
    final public Input<ParametricDistribution> distInput = new Input<>("distribution", "Distribution used to calculate likelihood", Validate.REQUIRED);
    // Optional placeholder if needed elsewhere; not used directly in this likelihood
    public final Input<RealParameter> uninfectiousRateInput = new Input<>(
            "uninfectiousRate", "Fixed uninfectious rate (per time unit), optional.", Validate.OPTIONAL);
    // Optional scaling of the mean (akin to a sampling rate in a given deme); defaults to 1.0 (no scaling)
    public final Input<RealParameter> scalingInput = new Input<>(
            "scaling", "Scaling factor applied to the prevalence-derived mean; must be > 0.", Validate.OPTIONAL);
    // Deme selection and filtering behavior
    public final Input<Integer> demeIndexInput = new Input<>(
            "demeIndex", "0-based deme index to match against CaseCountData traitIndices", Validate.REQUIRED);
    public final Input<Boolean> strictTraitFilteringInput = new Input<>(
            "strictTraitFiltering", "If true, warns when no observations match the specified deme; observations for other demes are ignored regardless.", Validate.OPTIONAL);
    
    protected CaseCountData caseCountData;
    protected ParametricDistribution dist;
    protected PrevalenceSkygrowth prevalence;
    protected int demeIndex;
    protected int[] filteredObservationIndices;
    protected boolean strictTraitFiltering = true;
    
    @Override
    public void initAndValidate() {
        caseCountData = caseCountInput.get();
        dist = distInput.get();
        prevalence = prevalenceSingleInput.get();
        Integer di = demeIndexInput.get();
        if (di == null) {
            throw new IllegalArgumentException("CaseCountLikelihood: 'demeIndex' must be provided and 0-based.");
        }
        demeIndex = di.intValue();
        Boolean stf = strictTraitFilteringInput.get();
        if (stf != null) {
            strictTraitFiltering = stf.booleanValue();
        }

        // Validate optional scaling parameter if present
        RealParameter scaleParam = scalingInput.get();
        if (scaleParam != null) {
            double s = scaleParam.getArrayValue();
            if (!(s > 0.0)) {
                throw new IllegalArgumentException("CaseCountLikelihood: 'scaling' must be > 0, got " + s);
            }
        }

        // Build filtered observation index list for this deme
        int n = caseCountData.getObservationCount();
        ArrayList<Integer> idxList = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            if (caseCountData.getTrait(i) == demeIndex) {
                idxList.add(i);
            }
        }
        if (idxList.isEmpty() && strictTraitFiltering) {
            System.err.println("Warning: CaseCountLikelihood found no observations for demeIndex=" + demeIndex + ". Likelihood will be neutral (0).");
        }
        filteredObservationIndices = new int[idxList.size()];
        for (int k = 0; k < idxList.size(); k++) filteredObservationIndices[k] = idxList.get(k);

        calculateLogP();
    }
    
    @Override
    public double calculateLogP() {
        logP = 0.0;
        // Neutral likelihood if there are no observations for this deme
        if (filteredObservationIndices == null || filteredObservationIndices.length == 0) {
            return logP;
        }

        // Calculate likelihood for each observation matching the selected deme
        for (int idx : filteredObservationIndices) {
            double t = caseCountData.getTime(idx);
            double caseCount = caseCountData.getCaseCount(idx);

            // Get log-prevalence at this time for the selected deme and convert to prevalence I
            double logI = prevalence.getPrevalenceTime(t);
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

            // Set the mean parameter for the distribution (based on prevalence I)
            if (dist instanceof GammaPoisson) {
                RealParameter meanParam = new RealParameter(new Double[]{scaledMean});
                ((GammaPoisson) dist).meanInput.setValue(meanParam, dist);
            } else {
                throw new IllegalArgumentException(
                        "CaseCountLikelihood currently supports only GammaPoisson as 'distribution'. Got: "
                                + dist.getClass().getName());
            }

            // Calculate log likelihood using the specified distribution
            double intervalLogP = dist.calcLogP(new RealParameter(new Double[]{caseCount}));
            if (Double.isNaN(intervalLogP)) {
                System.err.println("Warning: NaN likelihood for observation index " + idx + ": caseCount=" + caseCount + ", scaled mean=" + scaledMean);
                return Double.NEGATIVE_INFINITY;
            }

            logP += intervalLogP;
        }

        return logP;
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
