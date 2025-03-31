package mascot.distribution;

import org.apache.commons.math.distribution.PascalDistributionImpl;
import org.apache.commons.math.distribution.Distribution;
import org.apache.commons.math.distribution.PascalDistribution;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.parameter.RealParameter;

/**
 * Implementation of the Negative Binomial distribution (also known as Pascal distribution)
 * 
 * The probability mass function is given by
 * P(X = k) = C(k + r - 1, r - 1) * p^r * (1 - p)^k,
 * where r is the number of successes, p is the probability of success, and X is the total number of failures.
 * 
 * The mean and variance of X are
 * E(X) = (1 - p) * r / p, var(X) = (1 - p) * r / p^2.
 * Finally, the cumulative distribution function is given by
 * P(X <= k) = I(p, r, k + 1), where I is the regularized incomplete Beta function. 
 * 
 * Parameterized by mean and dispersion (alpha)
 * Where r = 1/alpha (number of successes)
 * and p = 1/(mean*alpha + 1) (probability of success)
 */
@Description("A negative binomial distribution parameterized by mean and dispersion.")
public class NegativeBinomialDistribution extends ParametricDistribution {
    final public Input<Function> meanInput = new Input<>("mean", "Mean parameter of the negative binomial distribution.");
    final public Input<Function> dispersionInput = new Input<>("dispersion", "Dispersion parameter (alpha) of the negative binomial distribution.");

    private PascalDistribution dist;

    // Must provide empty constructor for construction by XML
    public NegativeBinomialDistribution() {
        dist = new PascalDistributionImpl(1, 0.5); // Default values, will be updated in refresh()
    }

    public NegativeBinomialDistribution(RealParameter mean, RealParameter dispersion) {
        this();
        try {
            initByName("mean", mean, "dispersion", dispersion);
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException("Failed to initialize NegativeBinomialDistribution with mean and dispersion parameters.");
        }
    }

    @Override
    public void initAndValidate() {
        if (meanInput.get() != null && meanInput.get() instanceof RealParameter) {
            RealParameter meanParam = (RealParameter) meanInput.get();
            if (meanParam.getLower() == null) {
                meanParam.setLower(0.0);
            }
            if (meanParam.getUpper() == null) {
                meanParam.setUpper(Double.POSITIVE_INFINITY);
            }
        }

        if (dispersionInput.get() != null && dispersionInput.get() instanceof RealParameter) {
            RealParameter dispersionParam = (RealParameter) dispersionInput.get();
            if (dispersionParam.getLower() == null) {
                dispersionParam.setLower(0.0);
            }
            if (dispersionParam.getUpper() == null) {
                dispersionParam.setUpper(Double.POSITIVE_INFINITY);
            }
        }
        refresh();
    }

    /**
     * Make sure internal state is up to date
     */
    @SuppressWarnings("deprecation")
    void refresh() {
        double mean = meanInput.get() == null ? 1.0 : meanInput.get().getArrayValue();
        double dispersion = dispersionInput.get() == null ? 1.0 : dispersionInput.get().getArrayValue();
        
        if (mean <= 0) {
            throw new IllegalArgumentException("Mean parameter must be positive, got " + mean);
        }
        if (dispersion <= 0) {
            throw new IllegalArgumentException("Dispersion parameter must be positive, got " + dispersion);
        }
        
        // Convert mean and dispersion to r (number of successes) and p (probability of success)
        // The Pascal distribution is a special case of the Negative Binomial distribution where the number of successes parameter is an integer
        int r = (int) Math.round(1.0 / dispersion);
        if (r < 1) r = 1;
        double p = 1.0 / (mean * dispersion + 1.0);
        
        dist = new PascalDistributionImpl(r, p);
    }

    @Override
    public Distribution getDistribution() {
        refresh();
        return dist;
    }
    
    @Override
    public double getMeanWithoutOffset() {
        return meanInput.get().getArrayValue();
    }


}
