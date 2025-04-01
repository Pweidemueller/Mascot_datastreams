package mascot.distribution;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.parameter.RealParameter;
import mascot.parameterdynamics.NeDynamicsList;

import java.util.List;
import java.util.Random;

@Description("Calculates likelihood of observing case counts given Ne values using a specified distribution")
public class CaseCountLikelihood extends Distribution {
    public Input<NeDynamicsList> parametricFunctionInput = new Input<>(
            "NeDynamics", "input of the log effective population sizes", Validate.REQUIRED);  
    final public Input<CaseCountData> caseCountInput = new Input<>("caseCounts", "Case count data", Validate.REQUIRED);
    final public Input<ParametricDistribution> distInput = new Input<>("distribution", "Distribution used to calculate likelihood", Validate.REQUIRED);
    
    protected CaseCountData caseCountData;
    protected ParametricDistribution dist;
    protected NeDynamicsList parametricFunction;
    
    @Override
    public void initAndValidate() {
        caseCountData = caseCountInput.get();
        dist = distInput.get();
        parametricFunction = parametricFunctionInput.get();
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
            
            // Get Ne for this trait at this time
            double neValue = parametricFunction.get(trait).getNeTime(t);
            
            // Validate parameters
            if (neValue <= 0.0 || caseCount < 0.0) {
                System.err.println("Warning: Invalid parameters - Ne: " + neValue + ", caseCount: " + caseCount);
                return Double.NEGATIVE_INFINITY;
            }
            
            // Set the mean parameter for the distribution (based on Ne)
            if (dist instanceof NegativeBinomialDistribution) {
                RealParameter meanParam = new RealParameter(new Double[]{neValue});
                ((NegativeBinomialDistribution) dist).meanInput.setValue(meanParam, dist);
            }
            
            // Calculate log likelihood using the specified distribution
            double intervalLogP = dist.calcLogP(new RealParameter(new Double[]{caseCount}));
            if (Double.isNaN(intervalLogP)) {
                System.err.println("Warning: NaN likelihood for observation " + i + ": caseCount=" + caseCount + ", Ne=" + neValue);
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
