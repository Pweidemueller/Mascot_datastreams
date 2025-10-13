package mascotdatastreams.distribution;

import org.apache.commons.math.distribution.AbstractIntegerDistribution;
import org.apache.commons.math.distribution.Distribution;
import org.apache.commons.math.special.Beta;
import org.apache.commons.math.special.Gamma;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.parameter.RealParameter;

/**
 * Gamma-Poisson (Negative Binomial) distribution parameterised by mean (mu) and dispersion (alpha).
 *
 * Parameterization:
 *   r = 1 / alpha
 *   p = r / (r + mu)
 * so that E[X] = mu and Var[X] = mu + alpha * mu^2.
 *
 * PMF:
 *   P(X = k) = Γ(r + k) / (Γ(k+1) * Γ(r)) * p^r * (1 - p)^k
 * CDF:
 *   P(X <= k) = I(p; r, k + 1), regularized incomplete beta
 */
@Description("Gamma-Poisson (Negative Binomial) distribution parameterised by mean and dispersion (alpha).")
public class GammaPoisson extends ParametricDistribution implements CountDistributionWithMean {
    public final Input<Function> meanInput = new Input<>("mean", "Mean (mu) of the distribution.");
    public final Input<Function> dispersionInput = new Input<>("dispersion", "Dispersion (alpha) parameter.");

    private GammaPoissonDistributionImpl dist;

    // Empty constructor for XML
    public GammaPoisson() {
        dist = new GammaPoissonDistributionImpl(1.0, 0.5); // placeholder, refreshed in init/refresh
    }

    public GammaPoisson(RealParameter mean, RealParameter dispersion) {
        this();
        try {
            initByName("mean", mean, "dispersion", dispersion);
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException("Failed to initialize GammaPoisson with mean and dispersion parameters.");
        }
    }

    @Override
    public void initAndValidate() {
        if (meanInput.get() instanceof RealParameter) {
            RealParameter meanParam = (RealParameter) meanInput.get();
            if (meanParam.getLower() == null) {
                meanParam.setLower(0.0);
            }
            if (meanParam.getUpper() == null) {
                meanParam.setUpper(Double.POSITIVE_INFINITY);
            }
        }
        if (dispersionInput.get() instanceof RealParameter) {
            RealParameter dispParam = (RealParameter) dispersionInput.get();
            if (dispParam.getLower() == null) {
                dispParam.setLower(0.0);
            }
            if (dispParam.getUpper() == null) {
                dispParam.setUpper(Double.POSITIVE_INFINITY);
            }
        }
        refresh();
    }

    // Keep internal state up-to-date
    void refresh() {
        double mean = meanInput.get() == null ? 1.0 : meanInput.get().getArrayValue();
        double alpha = dispersionInput.get() == null ? 1.0 : dispersionInput.get().getArrayValue();
        if (mean <= 0) {
            throw new IllegalArgumentException("Mean must be > 0, got " + mean);
        }
        if (alpha <= 0) {
            throw new IllegalArgumentException("Dispersion (alpha) must be > 0, got " + alpha);
        }
        double r = 1.0 / alpha;
        double p = r / (r + mean);
        // Defensive clamp to avoid log(0)
        p = Math.min(1 - 1e-16, Math.max(1e-16, p));
        dist = new GammaPoissonDistributionImpl(r, p);
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

    // Public convenience methods requested: PMF and logPMF
    public double pmf(int x) {
        refresh();
        return dist.probability(x);
    }

    public double logPmf(int x) {
        refresh();
        return dist.logProbability(x);
    }

    @Override
    public double logPForMean(double observation, double mean) {
        // Stateless log PMF for a given mean; do not mutate Inputs
        if (observation < 0) {
            return Double.NEGATIVE_INFINITY;
        }
        // Dispersion (alpha) from input; default to 1.0 if absent, as in refresh()
        double alpha = dispersionInput.get() == null ? 1.0 : dispersionInput.get().getArrayValue();
        if (!(mean > 0.0) || !(alpha > 0.0)) {
            throw new IllegalArgumentException("GammaPoisson.logPForMean: mean and alpha must be > 0. Got mean=" + mean + ", alpha=" + alpha);
        }
        int x = (int) Math.round(observation);
        if (x < 0) {
            return Double.NEGATIVE_INFINITY;
        }
        double r = 1.0 / alpha;
        double p = r / (r + mean);
        p = Math.min(1 - 1e-16, Math.max(1e-16, p));
        // log PMF = ln Γ(r + x) - ln Γ(r) - ln Γ(x+1) + r ln p + x ln(1-p)
        double logGammaRK = Gamma.logGamma(r + x);
        double logGammaR = Gamma.logGamma(r);
        double logGammaK1 = Gamma.logGamma(x + 1.0);
        double logP = Math.log(p);
        double log1mP = Math.log(1.0 - p);
        return logGammaRK - logGammaR - logGammaK1 + r * logP + x * log1mP;
    }
}

/**
 * Inner integer distribution implementation.
 */
class GammaPoissonDistributionImpl extends AbstractIntegerDistribution {
    private static final long serialVersionUID = 1L;
	private final double r;        // shape
    private final double p;        // success probability
    private final double mean;     // derived mean
    private final double logGammaR;
    private final double logP;
    private final double log1mP;

    public GammaPoissonDistributionImpl(double r, double p) {
        if (!(r > 0.0)) throw new IllegalArgumentException("r must be > 0");
        if (!(p > 0.0 && p < 1.0)) throw new IllegalArgumentException("p must be in (0,1)");
        this.r = r;
        this.p = p;
        this.mean = r * (1.0 - p) / p;
        this.logGammaR = Gamma.logGamma(r);
        this.logP = Math.log(p);
        this.log1mP = Math.log(1.0 - p);
    }

    @Override
    public double probability(int x) {
        if (x < 0) return 0.0;
        return Math.exp(logProbability(x));
    }

    public double logProbability(int x) {
        if (x < 0) return Double.NEGATIVE_INFINITY;
        // log PMF = ln Γ(r + x) - ln Γ(r) - ln Γ(x+1) + r ln p + x ln(1-p)
        double logGammaRK = Gamma.logGamma(r + x);
        double logGammaK1 = Gamma.logGamma(x + 1.0);
        return logGammaRK - logGammaR - logGammaK1 + r * logP + x * log1mP;
    }

    @Override
    public double cumulativeProbability(int x) throws org.apache.commons.math.MathException {
        if (x < 0) return 0.0;
        // CDF: I(p; r, x+1)
        return Beta.regularizedBeta(p, r, x + 1);
    }

    @Override
    protected int getDomainLowerBound(double pQuantile) {
        return 0;
    }

    @Override
    protected int getDomainUpperBound(double pQuantile) {
        return Integer.MAX_VALUE;
    }

    public double getNumericalMean() {
        return mean;
    }

    public double getNumericalVariance() {
        return mean + (mean * mean) * (1.0 / r);
    }

    public double getSupportLowerBound() {
        return 0;
    }

    public double getSupportUpperBound() {
        return Double.POSITIVE_INFINITY;
    }

    public boolean isSupportConnected() {
        return true;
    }
}
