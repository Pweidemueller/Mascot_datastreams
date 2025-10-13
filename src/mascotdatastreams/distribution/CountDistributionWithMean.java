package mascotdatastreams.distribution;

/**
 * A small interface for count distributions that can compute a per-observation
 * log-likelihood given a mean without mutating BEAST Inputs during evaluation.
 */
public interface CountDistributionWithMean {
    /**
     * Compute log P(X = observation | mean, other distribution parameters)
     * without mutating any BEAST Inputs.
     */
    double logPForMean(double observation, double mean);
}
