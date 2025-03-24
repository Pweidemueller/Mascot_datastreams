package mascot.distribution;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import mascot.dynamics.Dynamics;

import java.util.List;
import java.util.Random;

import org.apache.commons.math.special.Gamma;

@Description("Calculates likelihood of observing case counts given Ne values using negative binomial distribution")
public class CaseCountLikelihood extends Distribution {
    final public Input<Dynamics> dynamicsInput = new Input<>("dynamics", "Input of Ne dynamics", Validate.REQUIRED);
    final public Input<CaseCountData> caseCountInput = new Input<>("caseCounts", "Case count data", Validate.REQUIRED);
    final public Input<RealParameter> dispersionInput = new Input<>("dispersion", "Dispersion parameter for negative binomial distribution", Validate.REQUIRED);
    
    protected CaseCountData caseCountData;
    protected double dispersion;
    protected Dynamics dynamics;
    protected double storedLogP;
    
    @Override
    public void initAndValidate() {
        caseCountData = caseCountInput.get();
        dispersion = dispersionInput.get().getValue();
        dynamics = dynamicsInput.get();
        calculateLogP();
    }
    
    @Override
    public double calculateLogP() {
        storedLogP = 0.0;
        
        // Calculate likelihood for each observation
        for (int i = 0; i < caseCountData.getObservationCount(); i++) {
            double time = caseCountData.getTime(i);
            int trait = caseCountData.getTrait(i);
            double caseCount = caseCountData.getCaseCount(i);
            
            // Find the correct interval for this time point
            int interval = 0;
            double currentTime = 0.0;
            while (interval < dynamics.getDimension() && currentTime + dynamics.getInterval(interval) < time) {
                currentTime += dynamics.getInterval(interval);
                interval++;
            } // i don't want Interval, should be some function that looks up Ne for a specific time: getNeTime(double t) (t is relative to most recent sample)
            
            // Get coalescent rate for this trait at this time
            // Note: coalescent rate is 1/Ne, so we need to invert it
            // TODO:
            // 1. when we are talking Ne do we have to do 1/coalrate -> Yes!
            // 2. given a date of a casecount how do I know "which" Ne segement sto request? given a trait how do I access the NE for that trait? -> time is releative to the most recent sample
            // 3. how do I put an operator on the dispersion parameter?
            // 4. don't use coalRates, but Ne directly!
            double[] coalRates = dynamics.getCoalescentRate(interval);
            if (coalRates == null || coalRates.length <= trait || coalRates[trait] <= 0.0) {
                System.err.println("Warning: Invalid coalescent rate at time " + time + ", trait " + trait + ", interval " + interval);
                return Double.NEGATIVE_INFINITY;
            }
            
            // Ne = 1/coalRate
            double neValue = 1.0 / coalRates[trait];
            
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
            
            storedLogP += intervalLogP;
        }
        
        logP = storedLogP;
        return storedLogP;
    }
    
    private double calculateNegBinomLogP(double x, double mean, double dispersion) {
        // Negative binomial log probability calculation
        // P(X = k) = Gamma(k + r)/(Gamma(r) * k!) * p^r * (1-p)^k
        // where r = dispersion, p = dispersion/(dispersion + mean)
        double p = dispersion / (dispersion + mean);
        double r = dispersion;
        
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
