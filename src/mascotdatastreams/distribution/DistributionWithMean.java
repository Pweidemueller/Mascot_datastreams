package mascotdatastreams.distribution;

/**
 * Interface for distributions parameterised by a mean, allowing stateless per-observation
 * log-likelihood evaluation without mutating BEAST Inputs.
 *
 * The {@code mean} parameter is always the mean of the distribution in its natural
 * parameterisation:
 * <ul>
 *   <li>{@link LogNormal}: {@code mean} is μ, the mean of the underlying normal distribution
 *       in log space (so the real-space median is exp(μ)).</li>
 *   <li>{@link GammaPoisson}: {@code mean} is the real-space count mean.</li>
 * </ul>
 */
public interface DistributionWithMean {
    /**
     * Compute log P(X = observation | mean, other distribution parameters)
     * without mutating any BEAST Inputs.
     */
    double logPForMean(double observation, double mean);
}
