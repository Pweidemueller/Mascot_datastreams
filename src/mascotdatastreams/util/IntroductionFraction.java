package mascotdatastreams.util;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Loggable;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import mascot.dynamics.RateShifts;
import mascotdatastreams.dynamics.Spline;


/**
 * Fraction of cases in a focal deme i that are direct introductions, integrated over the
 * window [recentTime, oldestTime] (backward time, i.e. time before the most recent sample):
 *
 *   P_i = int sum_{j != i} m^fw_{j->i} I_j(t) dt / int beta_i(t) I_i(t) dt
 *
 * The numerator is the cumulative import flux into deme i, the denominator the cumulative
 * total incidence of deme i (same beta * I integrand as SeroprevalenceLikelihood).
 * Without recentTime / oldestTime, the window defaults to the smallest / largest rate shift
 * (knot) of the focal spline.
 * Implements Function so it can be used as the "x" of a standard BEAST Prior.
 */
@Description("Cumulative fraction of incidence in a focal deme that is due to direct introductions "
        + "from other demes, integrated between recentTime and oldestTime (time before the most "
        + "recent sample; defaults to the smallest and largest rate shift of the focal spline). "
        + "Implements Function so that a prior can be placed on it.")
public class IntroductionFraction extends CalculationNode implements Function, Loggable {

    final public Input<Spline> splineInput = new Input<>("spline",
            "Prevalence spline of the focal deme i", Validate.REQUIRED);

    final public Input<List<Spline>> otherSplinesInput = new Input<>("otherSpline",
            "Prevalence splines of the source demes j != i", new ArrayList<>());

    final public Input<RealParameter> forwardMigrationInput = new Input<>("forwardMigration",
            "Forward migration rates m^fw_{j -> i}. Without forwardMigrationIndices, the dimension "
                    + "must equal the number of otherSpline entries (value k applies to otherSpline[k]).",
            Validate.REQUIRED);

    final public Input<IntegerParameter> forwardMigrationIndicesInput = new Input<>("forwardMigrationIndices",
            "Optional indices into forwardMigration, one per otherSpline entry (same order). Lets this "
                    + "object read from the shared flat migration vector of StructuredSkylinePrevalence, "
                    + "where the rate j -> i sits at index j*(n-1) + (i < j ? i : i-1).",
            Validate.OPTIONAL);

    final public Input<RealParameter> recentTimeInput = new Input<>("recentTime",
            "Most recent end of the integration window, as time before the most recent sample. "
                    + "Defaults to the smallest rate shift of the focal spline.",
            Validate.OPTIONAL);

    final public Input<RealParameter> oldestTimeInput = new Input<>("oldestTime",
            "Oldest end of the integration window, as time before the most recent sample. "
                    + "Defaults to the largest rate shift of the focal spline.",
            Validate.OPTIONAL);

    private Spline spline;
    private List<Spline> otherSplines;
    private RealParameter forwardMigration;
    private int[] forwardMigrationIndices;

    private double fraction;
    private boolean fractionKnown = false;

    @Override
    public void initAndValidate() {
        spline = splineInput.get();
        otherSplines = otherSplinesInput.get();
        forwardMigration = forwardMigrationInput.get();

        if (otherSplines.isEmpty()) {
            throw new IllegalArgumentException("IntroductionFraction: at least one otherSpline is required");
        }

        int n = otherSplines.size();
        forwardMigrationIndices = new int[n];
        IntegerParameter idxParam = forwardMigrationIndicesInput.get();
        if (idxParam == null) {
            if (forwardMigration.getDimension() != n) {
                throw new IllegalArgumentException("IntroductionFraction: forwardMigration dimension ("
                        + forwardMigration.getDimension() + ") must match the number of otherSpline entries ("
                        + n + "), or forwardMigrationIndices must be provided");
            }
            for (int k = 0; k < n; k++) forwardMigrationIndices[k] = k;
        } else {
            if (idxParam.getDimension() != n) {
                throw new IllegalArgumentException("IntroductionFraction: forwardMigrationIndices dimension ("
                        + idxParam.getDimension() + ") must match the number of otherSpline entries (" + n + ")");
            }
            for (int k = 0; k < n; k++) {
                int idx = idxParam.getValue(k);
                if (idx < 0 || idx >= forwardMigration.getDimension()) {
                    throw new IllegalArgumentException("IntroductionFraction: forwardMigrationIndices[" + k
                            + "] = " + idx + " is out of bounds for forwardMigration of dimension "
                            + forwardMigration.getDimension());
                }
                forwardMigrationIndices[k] = idx;
            }
        }

        if (!isWindowValid()) {
            RateShifts knots = spline.rateShiftsInput.get();
            throw new IllegalArgumentException("IntroductionFraction: window [recentTime, oldestTime] = ["
                    + getRecentTime() + ", " + getOldestTime() + "] must satisfy recentTime < oldestTime "
                    + "and lie within the rate shifts of the focal spline [" + knots.getValue(0) + ", "
                    + knots.getValue(knots.getDimension() - 1) + "]");
        }
        fractionKnown = false;
    }

    /**
     * True if recentTime < oldestTime and the window lies within the knots of the focal spline.
     * Outside the knots the spline is clamped to its boundary value, so the integral there would
     * not reflect the estimated dynamics.
     */
    private boolean isWindowValid() {
        RateShifts knots = spline.rateShiftsInput.get();
        double recentTime = getRecentTime();
        double oldestTime = getOldestTime();
        return recentTime < oldestTime
                && recentTime >= knots.getValue(0)
                && oldestTime <= knots.getValue(knots.getDimension() - 1);
    }

    /** Most recent end of the integration window (time before the most recent sample). */
    public double getRecentTime() {
        if (recentTimeInput.get() != null) return recentTimeInput.get().getArrayValue();
        return spline.rateShiftsInput.get().getValue(0);
    }

    /** Oldest end of the integration window (time before the most recent sample). */
    public double getOldestTime() {
        if (oldestTimeInput.get() != null) return oldestTimeInput.get().getArrayValue();
        RateShifts knots = spline.rateShiftsInput.get();
        return knots.getValue(knots.getDimension() - 1);
    }

    /** Import flux into the focal deme at time t: sum_j m^fw_{j->i} I_j(t). */
    public double getIntroductionIncidence(double t) {
        double sum = 0.0;
        for (int k = 0; k < otherSplines.size(); k++) {
            double m_ji = forwardMigration.getArrayValue(forwardMigrationIndices[k]);
            sum += m_ji * otherSplines.get(k).getPrevalence(t);
        }
        return sum;
    }

    /** Total incidence of the focal deme at time t: beta_i(t) I_i(t). */
    public double getTotalIncidence(double t) {
        return spline.getTransmissionRate(t) * spline.getPrevalence(t);
    }

    /**
     * Integrates import flux and total incidence with the trapezoid rule over the spline grid
     * points inside [recentTime, oldestTime], with the window ends as exact end points.
     * Times are backward (time before the most recent sample), so the loop runs from the most
     * recent end into the past; the direction does not affect the integrals.
     *
     * @return {cumulative introductions, cumulative total incidence}
     */
    public double[] getCumulativeIncidences() {
        double recentTime = getRecentTime();
        double oldestTime = getOldestTime();

        List<Double> times = new ArrayList<>();
        times.add(recentTime);
        for (int i = 0; i < spline.getGridPointCount(); i++) {
            double t = spline.getGridPointTime(i);
            if (t > recentTime && t < oldestTime) times.add(t);
        }
        times.add(oldestTime);

        double intro = 0.0;
        double total = 0.0;
        double prevT = times.get(0);
        double prevIntro = getIntroductionIncidence(prevT);
        double prevTotal = getTotalIncidence(prevT);
        for (int k = 1; k < times.size(); k++) {
            double t = times.get(k);
            double curIntro = getIntroductionIncidence(t);
            double curTotal = getTotalIncidence(t);
            double dt = t - prevT;
            intro += 0.5 * dt * (prevIntro + curIntro);
            total += 0.5 * dt * (prevTotal + curTotal);
            prevT = t;
            prevIntro = curIntro;
            prevTotal = curTotal;
        }
        return new double[] { intro, total };
    }

    private double calculateFraction() {
        if (!spline.update()) return Double.POSITIVE_INFINITY;
        // Only reachable if the window bounds are estimated and a proposal moves them out of range.
        if (!isWindowValid()) return Double.POSITIVE_INFINITY;
        double[] cum = getCumulativeIncidences();
        // Non-positive total incidence (possible with clipTransRate="false") has no meaningful
        // fraction; return a value outside [0, 1] so a bounded prior (e.g. Beta) rejects the state.
        if (!(cum[1] > 0.0)) return Double.POSITIVE_INFINITY;
        return cum[0] / cum[1];
    }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double getArrayValue() {
        if (!fractionKnown) {
            fraction = calculateFraction();
            fractionKnown = true;
        }
        return fraction;
    }

    @Override
    public double getArrayValue(int dim) {
        return getArrayValue();
    }

    @Override
    protected boolean requiresRecalculation() {
        // Splines always report dirty, so recompute whenever any upstream node changes.
        fractionKnown = false;
        return true;
    }

    @Override
    public void restore() {
        // The splines rebuild themselves after a restore, so recompute lazily from them.
        fractionKnown = false;
        super.restore();
    }

    @Override
    public void init(PrintStream out) {
        out.print(getID() == null ? "introductionFraction\t" : getID() + "\t");
    }

    @Override
    public void log(long sample, PrintStream out) {
        // Recompute rather than use the cache: when only logged (not part of the posterior),
        // requiresRecalculation() is never called and the cached value would never update.
        out.print(calculateFraction() + "\t");
    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }
}
