package mascot.distribution;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import mascot.parameterdynamics.NeDynamicsList;

import java.util.List;
import java.util.Random;

import org.apache.commons.math.special.Gamma;

@Description("Calculates likelihood of observing case counts given Ne values using negative binomial distribution")
public class CaseCountLikelihood extends Distribution {
    // final public Input<Dynamics> dynamicsInput = new Input<>("dynamics", "Input of Ne dynamics", Validate.REQUIRED);
    public Input<NeDynamicsList> parametricFunctionInput = new Input<>(
    		"NeDynamics", "input of the log effective population sizes", Validate.REQUIRED);  
    final public Input<CaseCountData> caseCountInput = new Input<>("caseCounts", "Case count data", Validate.REQUIRED);
    final public Input<RealParameter> dispersionInput = new Input<>("dispersion", "Dispersion parameter for negative binomial distribution", Validate.REQUIRED);
    
    protected CaseCountData caseCountData;
    protected double dispersion;
    // protected Dynamics dynamics;
    protected NeDynamicsList parametricFunction;
    
    @Override
    public void initAndValidate() {
        caseCountData = caseCountInput.get();
        dispersion = dispersionInput.get().getValue();
        // dynamics = dynamicsInput.get();
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
            if (neValue <= 0.0 || dispersion <= 0.0 || caseCount < 0.0) {
                System.err.println("Warning: Invalid parameters - Ne: " + neValue + ", dispersion: " + dispersion + ", caseCount: " + caseCount);
                return Double.NEGATIVE_INFINITY;
            }
            
            // Calculate negative binomial log likelihood
            double intervalLogP = calculateNegBinomLogP(caseCount, neValue, dispersion);
            if (Double.isNaN(intervalLogP)) {
                System.err.println("Warning: NaN likelihood for observation " + i + ": caseCount=" + caseCount + ", Ne=" + neValue + ", dispersion=" + dispersion);
                return Double.NEGATIVE_INFINITY;
            }
            
            logP += intervalLogP;
        }
        
        return logP;
    }
    
    private double calculateNegBinomLogP(double x, double mean, double dispersion) {
        // Negative binomial probability calculation
        // We parameterise the negative binomial distribution as
        // P(X=k) = Binom(k+r-1, k) * p^r * (1-p)^k
        // where
        // r = # successes until experiment stops
        // p = success probability
        // k = number of failures
        // Given dispersion alpha and mean, we have
        // r = 1/alpha
        // p = 1/(mean*alpha+1)
        double p = 1.0 / (mean * dispersion + 1.0);
        double r = 1.0 / dispersion;
    
        return Gamma.logGamma(x + r) - Gamma.logGamma(r) - Gamma.logGamma(x + 1) + 
               r * Math.log(p) + x * Math.log(1 - p);
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
