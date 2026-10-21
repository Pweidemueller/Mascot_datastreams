package mascotdatastreams.util;


import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.*;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import mascotdatastreams.dynamics.Spline;


/**
 * Hard constraint on the fraction of cases in a focal deme i that are direct introductions,
 * evaluated separately for every interval between two neighbouring grid points of the focal spline:
 *
 *   P_i[g_k, g_k+1] = int sum_{j != i} m^fw_{j->i} I_j(t) dt / int beta_i(t) I_i(t) dt
 *
 * (trapezoid rule on each interval, same integrands as IntroductionFraction).
 * Returns 0 if P_i <= maxFraction in every interval and negative infinity otherwise.
 */
@Description("returns 0 if the fraction of incidence in a focal deme that is due to direct introductions "
        + "is <= maxFraction in every grid interval of the focal spline, and negative infinity otherwise")
public class IntroductionFractionSmallerThan extends Distribution {
    final public Input<Spline> splineInput = new Input<>("spline",
            "Prevalence spline of the focal deme i", Validate.REQUIRED);

    final public Input<List<Spline>> otherSplinesInput = new Input<>("otherSpline",
            "Prevalence splines of the source demes j != i", new ArrayList<>());

    final public Input<RealParameter> forwardMigrationInput = new Input<>("forwardMigration",
            "Forward migration rates m^fw_{j -> i}. Without forwardMigrationIndices, the dimension "
                    + "must equal the number of otherSpline entries (value k applies to otherSpline[k]).",
            Validate.REQUIRED);

    final public Input<IntegerParameter> forwardMigrationIndicesInput = new Input<>("forwardMigrationIndices",
            "Optional indices into forwardMigration, one per otherSpline entry (same order).",
            Validate.OPTIONAL);

    final public Input<Double> maxFractionInput = new Input<>("maxFraction",
            "Maximum allowed introduction fraction in any grid interval", 0.5);

    private Spline spline;
    private List<Spline> otherSplines;
    private RealParameter forwardMigration;
    private int[] forwardMigrationIndices;


    @Override
    public void initAndValidate() {
        spline = splineInput.get();
        otherSplines = otherSplinesInput.get();
        forwardMigration = forwardMigrationInput.get();

        if (otherSplines.isEmpty()) {
            throw new IllegalArgumentException("IntroductionFractionSmallerThan: at least one otherSpline is required");
        }

        int n = otherSplines.size();
        forwardMigrationIndices = new int[n];
        IntegerParameter idxParam = forwardMigrationIndicesInput.get();
        if (idxParam == null) {
            if (forwardMigration.getDimension() != n) {
                throw new IllegalArgumentException("IntroductionFractionSmallerThan: forwardMigration dimension ("
                        + forwardMigration.getDimension() + ") must match the number of otherSpline entries ("
                        + n + "), or forwardMigrationIndices must be provided");
            }
            for (int k = 0; k < n; k++) forwardMigrationIndices[k] = k;
        } else {
            if (idxParam.getDimension() != n) {
                throw new IllegalArgumentException("IntroductionFractionSmallerThan: forwardMigrationIndices dimension ("
                        + idxParam.getDimension() + ") must match the number of otherSpline entries (" + n + ")");
            }
            for (int k = 0; k < n; k++) {
                int idx = idxParam.getValue(k);
                if (idx < 0 || idx >= forwardMigration.getDimension()) {
                    throw new IllegalArgumentException("IntroductionFractionSmallerThan: forwardMigrationIndices[" + k
                            + "] = " + idx + " is out of bounds for forwardMigration of dimension "
                            + forwardMigration.getDimension());
                }
                forwardMigrationIndices[k] = idx;
            }
        }

        calculateLogP();
    }

    /** Import flux into the focal deme at time t: sum_j m^fw_{j->i} I_j(t). */
    private double getIntroductionIncidence(double t) {
        double sum = 0.0;
        for (int k = 0; k < otherSplines.size(); k++) {
            double m_ji = forwardMigration.getArrayValue(forwardMigrationIndices[k]);
            sum += m_ji * otherSplines.get(k).getPrevalence(t);
        }
        return sum;
    }

    /** Total incidence of the focal deme at time t: beta_i(t) I_i(t). */
    private double getTotalIncidence(double t) {
        return spline.getTransmissionRate(t) * spline.getPrevalence(t);
    }

    @Override
    public double calculateLogP() {
        // Update the splines to ensure prevalences and transmission rates are calculated
        if (!spline.update()) {
            logP = Double.NEGATIVE_INFINITY;
            return Double.NEGATIVE_INFINITY;
        }
        for (Spline other : otherSplines) {
            if (!other.update()) {
                logP = Double.NEGATIVE_INFINITY;
                return Double.NEGATIVE_INFINITY;
            }
        }

        double maxFraction = maxFractionInput.get();

        // Integrate each grid interval with the trapezoid rule and check its introduction fraction
        // this is going backwards in time, but since we just care about the number of introductions and new cases in an interval between two neighboring grid points, the direction doesn't matter
        int gridPointCount = spline.getGridPointCount();
        double prevT = spline.getGridPointTime(0);
        double prevIntro = getIntroductionIncidence(prevT);
        double prevTotal = getTotalIncidence(prevT);
        for (int i = 1; i < gridPointCount; i++) {
            double t = spline.getGridPointTime(i);
            double curIntro = getIntroductionIncidence(t);
            double curTotal = getTotalIncidence(t);
            double dt = t - prevT;
            double intro = 0.5 * dt * (prevIntro + curIntro);
            double total = 0.5 * dt * (prevTotal + curTotal);
            // Non-positive total incidence (possible with clipTransRate="false") has no
            // meaningful fraction, so it is rejected as well
            if (!(total > 0.0) || intro / total > maxFraction) {
                logP = Double.NEGATIVE_INFINITY;
                return Double.NEGATIVE_INFINITY;
            }
            prevT = t;
            prevIntro = curIntro;
            prevTotal = curTotal;
        }

        logP = 0.0;
        return 0.0;
    }


    @Override
    public void sample(State state, Random random) {
		// do nothing


    }

    @Override
    public List<String> getConditions() {
        List<String> conditions = new ArrayList<>();
        return conditions;
    }

    @Override
    public List<String> getArguments() {
        List<String> arguments = new ArrayList<>();

        return arguments;
    }
}
