package mascotdatastreams.distribution;

import org.apache.commons.math.distribution.AbstractIntegerDistribution;
import org.apache.commons.math.distribution.Distribution;
import org.apache.commons.math.special.Gamma;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.distribution.ParametricDistribution;

/**
 * Binomial distribution parameterised by number of trials (n) and success probability (p).
 *
 * PMF: P(X = k) = C(n, k) p^k (1 - p)^(n - k)
 */
@Description("Binomial distribution parameterised by trials (n) and probability (p).")
public class Binomial extends ParametricDistribution implements BinomialDistributionWithProb {
    private static final double PROB_EPS = 1e-16;

    public final Input<Function> nInput = new Input<>("n", "Number of trials (non-negative integer).");
    public final Input<Function> pInput = new Input<>("p", "Success probability in (0,1).");

    private BinomialDistributionImpl dist;

    public Binomial() {
        dist = new BinomialDistributionImpl(1, 0.5);
    }

    @Override
    public void initAndValidate() {
        refresh();
    }

    void refresh() {
        int n = 1;
        double p = 0.5;
        if (nInput.get() != null) {
            n = (int) Math.round(nInput.get().getArrayValue());
        }
        if (pInput.get() != null) {
            p = pInput.get().getArrayValue();
        }
        if (n < 0) {
            throw new IllegalArgumentException("Binomial: n must be >= 0, got " + n);
        }
        // clamp p to avoid log(0)
        p = Math.min(1.0 - PROB_EPS, Math.max(PROB_EPS, p));
        dist = new BinomialDistributionImpl(n, p);
    }

    @Override
    public Distribution getDistribution() {
        refresh();
        return dist;
    }

    public double pmf(int x) {
        refresh();
        return dist.probability(x);
    }

    public double logPmf(int x) {
        refresh();
        return dist.logProbability(x);
    }

    /**
     * Stateless log PMF for convenience.
     */
    public static double logPMF(int x, int n, double p) {
        return BinomialDistributionImpl.staticLogProbability(x, n, p);
    }

    /**
     * Instance method variant so callers can use an injected Binomial distribution
     * without mutating its inputs. Delegates to a stateless implementation.
     */
    @Override
    public double logPMFForParams(int x, int n, double p) {
        return BinomialDistributionImpl.staticLogProbability(x, n, p);
    }
}

class BinomialDistributionImpl extends AbstractIntegerDistribution {
    private static final long serialVersionUID = 1L;
    private static final double PROB_EPS = 1e-16;

    private final int n;
    private final double p;
    private final double logP;
    private final double log1mP;

    public BinomialDistributionImpl(int n, double p) {
        if (n < 0) throw new IllegalArgumentException("n must be >= 0");
        if (!(p > 0.0 && p < 1.0)) throw new IllegalArgumentException("p must be in (0,1)");
        this.n = n;
        this.p = p;
        this.logP = Math.log(p);
        this.log1mP = Math.log(1.0 - p);
    }

    @Override
    public double probability(int x) {
        if (x < 0 || x > n) return 0.0;
        return Math.exp(logProbability(x));
    }

    public double logProbability(int x) {
        return staticLogProbability(x, n, p);
    }

    static double staticLogProbability(int x, int n, double p) {
        if (n < 0) return Double.NEGATIVE_INFINITY;
        if (x < 0 || x > n) return Double.NEGATIVE_INFINITY;
        // clamp p for numerical stability
        double pp = Math.min(1.0 - PROB_EPS, Math.max(PROB_EPS, p));
        // log C(n,k) + k log p + (n-k) log (1-p)
        double logComb = Gamma.logGamma(n + 1.0) - Gamma.logGamma(x + 1.0) - Gamma.logGamma(n - x + 1.0);
        return logComb + x * Math.log(pp) + (n - x) * Math.log(1.0 - pp);
    }

    @Override
    public double cumulativeProbability(int x) throws org.apache.commons.math.MathException {
        // Simple summation; not optimized as we only need log-probability in this project
        if (x < 0) return 0.0;
        if (x >= n) return 1.0;
        double sum = 0.0;
        for (int k = 0; k <= x; k++) {
            sum += Math.exp(staticLogProbability(k, n, p));
        }
        return Math.min(1.0, Math.max(0.0, sum));
    }

    @Override
    protected int getDomainLowerBound(double pQuantile) {
        return 0;
    }

    @Override
    protected int getDomainUpperBound(double pQuantile) {
        return n;
    }
}


