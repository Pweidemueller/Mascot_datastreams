package mascotdatastreams.distribution;

/**
 * A small interface for distributions (both discrete and continuous) that can compute
 * a per-observation log-likelihood given a mean without mutating BEAST Inputs during evaluation.
 * 
 * This interface is used by both count distributions (e.g., GammaPoisson) and continuous
 * distributions (e.g., LogNormal) that can be parameterized by their mean.
 */
public interface DistributionWithMean {
    /**
     * Compute log P(X = observation | mean, other distribution parameters)
     * without mutating any BEAST Inputs.
     * 
     * For discrete distributions, this is log P(X = observation).
     * For continuous distributions, this is log f(observation) where f is the PDF.
     */
    double logPForMean(double observation, double mean);
}
