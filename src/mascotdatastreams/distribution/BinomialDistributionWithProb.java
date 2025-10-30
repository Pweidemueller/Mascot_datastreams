package mascotdatastreams.distribution;

/**
 * A small interface for count distributions that can compute a per-observation
 * log-likelihood given a mean without mutating BEAST Inputs during evaluation.
 */
public interface BinomialDistributionWithProb {
    /**
     * Compute log P(X = x | n, p)
     * without mutating any BEAST Inputs.
     */
    double logPMFForParams(int x, int n, double p);
}
