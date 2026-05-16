package mascotdatastreams.dynamics;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import mascot.parameterdynamics.NeDynamics;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Maps log-prevalence dynamics to an Ne(t) process using spline interpolation.
 * 
 * This class takes a Spline for log-prevalence interpolation and maps
 * the prevalence to coalescent effective population size Ne(t) using the transmission rate
 * formula: Ne(t) = I(t) / (c * transmission_rate(t)), where transmission_rate(t) 
 * is the derivative of log-prevalence plus the uninfectious rate.
 * 
 * Key differences from Skygrowth:
 * - Uses spline interpolation instead of piecewise-exponential
 * - Takes a pre-computed Spline as input rather than individual parameters
 * - Provides analytical derivatives via spline
 * - Uses precomputed grid points for efficient lookup
 */
@Description("Maps log-prevalence to Ne(t) using spline interpolation between rate shift points")
public class SplinePrevalenceToNe extends NeDynamics implements Loggable {
    
    // Required inputs
    final public Input<Spline> splineInput = new Input<>("spline",
            "Spline for log-prevalence interpolation", 
            Input.Validate.REQUIRED);
    
    // Optional inputs
    final public Input<RealParameter> coalescentScaleInput = new Input<>("coalescentScale",
            "Coalescent scaling constant c in Ne = I / (c * transmission_rate)", 
            Input.Validate.OPTIONAL);

    final public Input<RealParameter> NeScalerInput = new Input<>("NeScaler",
            "Additional scaling factor applied to Ne(t): Ne = NeScaler * I / (c * transmission_rate). Default: 1.0",
            Input.Validate.OPTIONAL);

    final public Input<List<Spline>> otherSplinesInput = new Input<>("otherSpline",
            "Splines for other demes (j != i). If provided together with incomingForwardMigration, "
                    + "the migration contribution to the transmission rate is subtracted to obtain "
                    + "the local transmission rate, which is then used in getNeTime.",
            new ArrayList<>());

    final public Input<RealParameter> incomingForwardMigrationInput = new Input<>("incomingForwardMigration",
            "Forward migration rates m^fw_{j -> i} into this deme. By default, the parameter's dimension "
                    + "must equal the number of otherSpline entries (value k applies to otherSpline[k]). "
                    + "If forwardMigrationIndices is also provided, this parameter can be a larger shared "
                    + "vector (e.g. the same flat migration matrix used by StructuredSkylinePrevalence) "
                    + "and only the selected indices are read. Required when otherSpline is non-empty.",
            Input.Validate.OPTIONAL);

    final public Input<IntegerParameter> forwardMigrationIndicesInput = new Input<>("forwardMigrationIndices",
            "Optional indices into incomingForwardMigration, one per otherSpline entry (same order). "
                    + "Lets multiple SplinePrevalenceToNe instances share a single flat migration parameter "
                    + "while each picks out the rates incoming into its own deme.",
            Input.Validate.OPTIONAL);

    // Member variables
    private Spline spline;
    private RealParameter coalescentScale;
    private RealParameter neScaler;
    private List<Spline> otherSplines;
    private RealParameter incomingForwardMigration;
    private int[] forwardMigrationIndices;
    private boolean hasMigration;
    
    boolean NesKnown = false;
    boolean returnNaN = false;
    
    
    @Override
    public void initAndValidate() {
        spline = splineInput.get();
        coalescentScale = coalescentScaleInput.get();
        neScaler = NeScalerInput.get();
        if (neScaler == null) {
            neScaler = new RealParameter(new Double[] { 1.0 });
        }
        
        // Set time flag for NeDynamics
        isTime = true;
        
        // Validate spline is not null
        if (spline == null) {
            throw new IllegalArgumentException("spline input is required");
        }

        otherSplines = otherSplinesInput.get();
        incomingForwardMigration = incomingForwardMigrationInput.get();
        IntegerParameter idxParam = forwardMigrationIndicesInput.get();
        hasMigration = otherSplines != null && !otherSplines.isEmpty();
        if (hasMigration) {
            if (incomingForwardMigration == null) {
                throw new IllegalArgumentException(
                        "incomingForwardMigration is required when otherSpline is provided");
            }
            int n = otherSplines.size();
            forwardMigrationIndices = new int[n];
            if (idxParam == null) {
                if (incomingForwardMigration.getDimension() != n) {
                    throw new IllegalArgumentException(
                            "incomingForwardMigration dimension (" + incomingForwardMigration.getDimension()
                                    + ") must match the number of otherSpline entries (" + n + "); "
                                    + "alternatively, provide forwardMigrationIndices to pick elements out of "
                                    + "a larger shared parameter.");
                }
                for (int k = 0; k < n; k++) forwardMigrationIndices[k] = k;
            } else {
                if (idxParam.getDimension() != n) {
                    throw new IllegalArgumentException(
                            "forwardMigrationIndices dimension (" + idxParam.getDimension()
                                    + ") must match the number of otherSpline entries (" + n + ")");
                }
                int dim = incomingForwardMigration.getDimension();
                for (int k = 0; k < n; k++) {
                    int idx = idxParam.getValue(k);
                    if (idx < 0 || idx >= dim) {
                        throw new IllegalArgumentException(
                                "forwardMigrationIndices[" + k + "] = " + idx
                                        + " is out of bounds for incomingForwardMigration of dimension " + dim);
                    }
                    forwardMigrationIndices[k] = idx;
                }
            }
        }
    }
    
    public List<String> getParameterIds() {
        return null;
    }
    
    /**
     * Gets prevalence at time t using precomputed grid points for efficiency.
     * Returns the prevalence value at the closest grid point.
     * 
     * @param t time (backward from present)
     * @return prevalence at time t
     */
    public double getPrevalenceTime(double t) {
        return spline.getPrevalence(t);
    }

    /**
     * Migration contribution to the transmission rate at time t:
     *   beta^migration_i(t) = sum_{j != i} m^fw_{j -> i} * I_j(t) / I_i(t).
     *
     * Returns 0 if no otherSpline/incomingForwardMigration inputs were provided.
     */
    public double getMigrationTransmissionRate(double t) {
        if (!hasMigration) return 0.0;
        double I_i = spline.getPrevalence(t);
        double sum = 0.0;
        for (int k = 0; k < otherSplines.size(); k++) {
            double m_ji = incomingForwardMigration.getArrayValue(forwardMigrationIndices[k]);
            double I_j = otherSplines.get(k).getPrevalence(t);
            sum += m_ji * I_j / I_i;
        }
        return sum;
    }

    /**
     * Local transmission rate at time t:
     *   beta^local_i(t) = beta_i(t) - beta^migration_i(t).
     *
     * Falls back to the full spline transmission rate when no migration inputs are set.
     */
    public double getLocalTransmissionRate(double t) {
        double beta = spline.getTransmissionRate(t);
        if (!hasMigration) return beta;
        return beta - getMigrationTransmissionRate(t);
    }

    // First non-positive transmission-rate event is reported to stderr; subsequent ones are
    // counted silently so MCMC runs that frequently visit such states don't flood the log.
    private static boolean nonPositiveRateWarned = false;
    private static long nonPositiveRateCount = 0;

    /** Number of times getNeTime has seen a non-positive effective transmission rate. */
    public static long getNonPositiveRateCount() {
        return nonPositiveRateCount;
    }

    @Override
    public double getNeTime(double t) {
        // Get prevalence using precomputed grid points for efficiency
        double I_t = getPrevalenceTime(t);

        // Use the local (intra-deme) transmission rate when migration info is provided,
        // otherwise fall back to the aggregate spline transmission rate.
        double effectiveRate = hasMigration ? getLocalTransmissionRate(t) : spline.getTransmissionRate(t);

        if (effectiveRate <= 0.0) {
            nonPositiveRateCount++;
            if (!nonPositiveRateWarned) {
                nonPositiveRateWarned = true;
                System.err.println("Warning: non-positive "
                        + (hasMigration ? "local " : "")
                        + "transmission rate at time " + t + ": " + effectiveRate
                        + ". Further occurrences will be counted silently; use "
                        + "SplinePrevalenceToNe.getNonPositiveRateCount() to read the total.");
            }
        }

        // Get coalescent scaling constant
        double c = coalescentScale.getArrayValue();
        double scaler = neScaler.getArrayValue();

        // Compute Ne: Ne = I / (c * effectiveRate)
        return scaler * I_t / (c * effectiveRate);
    }
    
    // TODO: change to more salient recalculation criterion, e.g. just recalculate the splines and then return super.requiresRecalculation()
    @Override
    public boolean requiresRecalculation() {
        return true;  // Always recalculate when parameters change
    }
    
    // TODO do more salient isDirty return
    @Override
    public boolean isDirty() {        
        
        return true;
    }
    
    @Override
    public void store() {
        super.store();
    }

    @Override
    public void restore() {
        super.restore();
    }

    @Override
    public void init(PrintStream printStream) {
        for (int i = 0; i < spline.getGridPointCount(); i+=10) {
            printStream.print("I_" + i + "\t");
        }
        for (int i = 0; i < spline.getGridPointCount(); i+=10) {
            printStream.print("Ne_" + i + "\t");
        }
        for (int i = 0; i < spline.getGridPointCount(); i+=10) {
            printStream.print("transmissionRate_" + i + "\t");
        }
        for (int i = 0; i < spline.getGridPointCount(); i+=10) {
            printStream.print("localTransmissionRate_" + i + "\t");
        }
    }

    @Override
    public void log(long l, PrintStream printStream) {

        for (int i = 0; i < spline.getGridPointCount(); i+=10) {
            double t = spline.getGridPointTime(i);
            double prevalence = getPrevalenceTime(t);
            printStream.print(Math.log(prevalence) + "\t");
        }
        for (int i = 0; i < spline.getGridPointCount(); i+=10) {
            double t = spline.getGridPointTime(i);
            double Ne = getNeTime(t);
            printStream.print(Math.log(Ne) + "\t");
        }
        for (int i = 0; i < spline.getGridPointCount(); i+=10) {
            double t = spline.getGridPointTime(i);
            double transmissionRate = spline.getTransmissionRate(t);
            printStream.print(transmissionRate + "\t");
        }
        for (int i = 0; i < spline.getGridPointCount(); i+=10) {
            double t = spline.getGridPointTime(i);
            printStream.print(getLocalTransmissionRate(t) + "\t");
        }
    }

    @Override
    public void close(PrintStream printStream) {
        // Nothing to close
    }
}