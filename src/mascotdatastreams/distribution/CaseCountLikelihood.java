package mascotdatastreams.distribution;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.parameter.RealParameter;
import mascot.parameterdynamics.NeDynamicsList;
import mascotdatastreams.dynamics.PrevalenceDynamicsList;

import java.util.List;
import java.util.Random;

@Description("Calculates likelihood of observing case counts given prevalence values using a specified distribution")
public class CaseCountLikelihood extends Distribution {
    // Deprecated legacy input to provide a clear migration error if still used
    public final Input<NeDynamicsList> deprecatedNeInput = new Input<>(
            "NeDynamics", "Deprecated: supply prevalence dynamics instead of NeDynamics.", Validate.OPTIONAL);

    // New prevalence input
    public final Input<PrevalenceDynamicsList> prevalenceInput = new Input<>(
            "prevalence", "Input list of log-prevalence dynamics per deme", Validate.REQUIRED);

    final public Input<CaseCountData> caseCountInput = new Input<>("caseCounts", "Case count data", Validate.REQUIRED);
    final public Input<ParametricDistribution> distInput = new Input<>("distribution", "Distribution used to calculate likelihood", Validate.REQUIRED);
    // Optional placeholder if needed elsewhere; not used directly in this likelihood
    public final Input<RealParameter> uninfectiousRateInput = new Input<>(
            "uninfectiousRate", "Fixed uninfectious rate (per time unit), optional.", Validate.OPTIONAL);
    
    protected CaseCountData caseCountData;
    protected ParametricDistribution dist;
    protected PrevalenceDynamicsList prevalence;
    
    @Override
    public void initAndValidate() {
        caseCountData = caseCountInput.get();
        dist = distInput.get();
        if (deprecatedNeInput.get() != null) {
            throw new IllegalArgumentException("CaseCountLikelihood: The 'NeDynamics' input is deprecated. Please provide 'prevalence' (log-prevalence dynamics) instead.");
        }
        prevalence = prevalenceInput.get();
        calculateLogP();
    }
    
    @Override
    public double calculateLogP() {
        logP = 0.0;
        
        // Calculate likelihood for each observation
        for (int i = 0; i < caseCountData.getObservationCount(); i++) {
            double t = caseCountData.getTime(i);
            int trait = caseCountData.getTrait(i);
            double caseCount = caseCountData.getCaseCount(i);
            
            // Get log-prevalence for this trait at this time and convert to prevalence I
            double logI = prevalence.get(trait).getPrevalenceTime(t);
            double meanI = Math.exp(logI);
            
            // Validate parameters
            if (meanI <= 0.0 || caseCount < 0.0) {
                System.err.println("Warning: Invalid parameters - prevalence I: " + meanI + ", caseCount: " + caseCount);
                return Double.NEGATIVE_INFINITY;
            }
            
            // Set the mean parameter for the distribution (based on prevalence I)
            if (dist instanceof GammaPoisson) {
                RealParameter meanParam = new RealParameter(new Double[]{meanI});
                ((GammaPoisson) dist).meanInput.setValue(meanParam, dist);
            } else if (dist instanceof NegativeBinomialDistribution) { // backward compatibility if present
                RealParameter meanParam = new RealParameter(new Double[]{meanI});
                ((NegativeBinomialDistribution) dist).meanInput.setValue(meanParam, dist);
            }
            
            // Calculate log likelihood using the specified distribution
            double intervalLogP = dist.calcLogP(new RealParameter(new Double[]{caseCount}));
            if (Double.isNaN(intervalLogP)) {
                System.err.println("Warning: NaN likelihood for observation " + i + ": caseCount=" + caseCount + ", prevalence I=" + meanI);
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
