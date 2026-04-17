package mascotdatastreams.distribution;

import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.Distribution;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.parameter.RealParameter;

/**
 * Log-normal distribution parameterised by μ (mean of log(X)) and σ (sd of log(X)).
 *
 * If X ~ LogNormal(μ, σ²) then log(X) ~ Normal(μ, σ²).
 *   median(X) = exp(μ),   E[X] = exp(μ + σ²/2)
 *
 * PDF:
 *   f(x) = 1/(x σ √(2π)) · exp(-½ ((log x − μ)² / σ²))
 */
@Description("Log-normal distribution parameterised by mean (on log scale) and standard deviation (on log scale).")
public class LogNormal extends ParametricDistribution implements DistributionWithMean {
    public final Input<Function> meanInput = new Input<>("mean", "Mean on log scale, i.e. μ = E[log(X)].");
    public final Input<Function> sdInput = new Input<>("sd", "Standard deviation on log scale, σ = sd(log(X)).");

    private LogNormalDistributionImpl dist;

    // Empty constructor for XML
    public LogNormal() {
        dist = new LogNormalDistributionImpl(0.0, 0.5); // placeholder, refreshed in init/refresh
    }

    public LogNormal(RealParameter mean, RealParameter sd) {
        this();
        try {
            initByName("mean", mean, "sd", sd);
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException("Failed to initialize LogNormal with mean and sd parameters.");
        }
    }

    @Override
    public void initAndValidate() {
        if (sdInput.get() instanceof RealParameter) {
            RealParameter sdParam = (RealParameter) sdInput.get();
            if (sdParam.getLower() == null) {
                sdParam.setLower(0.0);
            }
            if (sdParam.getUpper() == null) {
                sdParam.setUpper(Double.POSITIVE_INFINITY);
            }
        }
        refresh();
    }

    void refresh() {
        double mu = meanInput.get() == null ? 0.0 : meanInput.get().getArrayValue();
        double sd = sdInput.get() == null ? 1.0 : sdInput.get().getArrayValue();
        if (sd <= 0) {
            throw new IllegalArgumentException("Standard deviation must be > 0, got " + sd);
        }
        dist = new LogNormalDistributionImpl(mu, sd);
    }

    @Override
    public Distribution getDistribution() {
        refresh();
        return dist;
    }

    @Override
    public double getMeanWithoutOffset() {
        Function meanFn = meanInput.get();
        if (meanFn == null) {
            return 0.0; // consistent with refresh() default (μ on log scale)
        }
        return meanFn.getArrayValue();
    }

    // Public convenience methods
    public double pdf(double x) {
        refresh();
        return dist.density(x);
    }

    public double logPdf(double x) {
        refresh();
        return dist.logDensity(x);
    }

    @Override
    public double logPForMean(double observation, double mean) {
        // mean is μ, the mean of the underlying normal in log space (= log of the real-space median).
        // No σ-dependent conversion: the caller is responsible for computing μ = log(median) upstream.
        if (observation <= 0.0) {
            return Double.NEGATIVE_INFINITY;
        }
        double sd = sdInput.get() == null ? 1.0 : sdInput.get().getArrayValue();
        if (!(sd > 0.0)) {
            throw new IllegalArgumentException("LogNormal.logPForMean: sd must be > 0. Got sd=" + sd);
        }
        return new LogNormalDistributionImpl(mean, sd).logDensity(observation);
    }
}

/**
 * Inner continuous distribution implementation.
 */
class LogNormalDistributionImpl implements ContinuousDistribution {
    private static final long serialVersionUID = 1L;
    private static final double LOG_SQRT_2PI = Math.log(Math.sqrt(2.0 * Math.PI));

    private final double mu;        // mean of log(X)
    private final double sigma;     // standard deviation of log(X)
    private final double logSigma;   // log(sigma) for efficiency

    public LogNormalDistributionImpl(double mu, double sigma) {
        if (!(sigma > 0.0)) throw new IllegalArgumentException("sigma must be > 0");
        this.mu = mu;
        this.sigma = sigma;
        this.logSigma = Math.log(sigma);
    }

    @Override
    public double density(double x) {
        if (x <= 0.0) return 0.0;
        return Math.exp(logDensity(x));
    }

    public double logDensity(double x) {
        if (x <= 0.0) return Double.NEGATIVE_INFINITY;
        // log f(x) = -log(x) - log(σ) - log(2π) - 0.5*((log(x) - μ)² / σ²)
        double logX = Math.log(x);
        double diff = logX - mu;
        return -logX - logSigma - LOG_SQRT_2PI - 0.5 * (diff * diff) / (sigma * sigma);
    }

    @Override
    public double cumulativeProbability(double x) throws org.apache.commons.math.MathException {
        if (x <= 0.0) return 0.0;
        // CDF: Φ((log(x) - μ) / σ) where Φ is the standard normal CDF
        double z = (Math.log(x) - mu) / sigma;
        return 0.5 * (1.0 + org.apache.commons.math.special.Erf.erf(z / Math.sqrt(2.0)));
    }

    @Override
    public double cumulativeProbability(double x0, double x1) throws org.apache.commons.math.MathException {
        return cumulativeProbability(x1) - cumulativeProbability(x0);
    }

    @Override
    public double inverseCumulativeProbability(double p) throws org.apache.commons.math.MathException {
        // TODO: This method has not been validated yet. The implementation below is kept for future investigation.
        throw new UnsupportedOperationException(
            "LogNormal.inverseCumulativeProbability() is not yet implemented/validated. " +
            "The implementation using approximateNormalQuantile() exists but needs verification. " +
            "This method is not currently used by WastewaterLikelihood, which only requires logPForMean()."
        );
        
        // Implementation kept for later investigation:
        // Inverse CDF: exp(μ + σ * Φ⁻¹(p))
        // Using approximation for standard normal quantile
        // if (p <= 0.0) return 0.0;
        // if (p >= 1.0) return Double.POSITIVE_INFINITY;
        // Simple approximation: z ≈ sqrt(2) * erfinv(2p - 1)
        // For now, use a simple numerical approximation
        // double z = approximateNormalQuantile(p);
        // return Math.exp(mu + sigma * z);
    }
    
    private double approximateNormalQuantile(double p) {
        // TODO: This method has not been validated yet. Implementation kept for future investigation.
        // Beasley-Springer-Moro algorithm approximation for normal quantile
        // Simplified version for p in (0,1)
        if (p < 0.5) {
            double t = Math.sqrt(-2.0 * Math.log(p));
            return -(t - (2.515517 + 0.802853 * t + 0.010328 * t * t) / 
                    (1.0 + 1.432788 * t + 0.189269 * t * t + 0.001308 * t * t * t));
        } else {
            double t = Math.sqrt(-2.0 * Math.log(1.0 - p));
            return (t - (2.515517 + 0.802853 * t + 0.010328 * t * t) / 
                   (1.0 + 1.432788 * t + 0.189269 * t * t + 0.001308 * t * t * t));
        }
    }

    public double getNumericalMean() {
        // E[X] = exp(μ + σ²/2)
        return Math.exp(mu + 0.5 * sigma * sigma);
    }

    public double getNumericalVariance() {
        // Var[X] = (exp(σ²) - 1) * exp(2μ + σ²)
        double SigmaSq = sigma * sigma;
        return (Math.exp(SigmaSq) - 1.0) * Math.exp(2.0 * mu + SigmaSq);
    }

    public double getSupportLowerBound() {
        return 0.0;
    }

    public double getSupportUpperBound() {
        return Double.POSITIVE_INFINITY;
    }

    public boolean isSupportConnected() {
        return true;
    }
}
