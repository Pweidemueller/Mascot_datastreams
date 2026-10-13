package mascotdatastreams.util;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import org.junit.jupiter.api.Test;

import beast.base.inference.distribution.Beta;
import beast.base.inference.distribution.Prior;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import mascot.dynamics.RateShifts;
import mascotdatastreams.dynamics.Spline;

public class IntroductionFractionTest {

    private static final String KNOTS = "0.0 0.2 0.4 0.6 0.8 1.0";
    private static final double GAMMA = 10.0;

    private static RateShifts rateShifts(String values) {
        RateShifts rs = new RateShifts();
        rs.initByName("value", values);
        return rs;
    }

    private static String fineGrid(double end, int n) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i <= n; i++) {
            sb.append(end * i / n).append(" ");
        }
        return sb.toString().trim();
    }

    private static Spline spline(Double[] logI, RateShifts knots, RateShifts grid) {
        Spline s = new Spline();
        s.initByName(
                "logInfected", new RealParameter(logI),
                "rateShifts", knots,
                "gridRateShifts", grid,
                "uninfectiousRate", new RealParameter(new Double[] { GAMMA }),
                "clipTransRate", false);
        return s;
    }

    private static Double[] linearLogI(double logI0, double slope, RateShifts knots) {
        Double[] v = new Double[knots.getDimension()];
        for (int k = 0; k < v.length; k++) {
            v[k] = logI0 + slope * knots.getValue(k);
        }
        return v;
    }

    @Test
    public void testConstantPrevalence() {
        // Constant I_i and I_j: beta_i = gamma, so P_i = m * I_j / (gamma * I_i).
        RateShifts knots = rateShifts(KNOTS);
        RateShifts grid = rateShifts(fineGrid(1.0, 100));
        double Ii = 100.0;
        double Ij = 50.0;
        double m = 2.0;

        IntroductionFraction f = new IntroductionFraction();
        f.initByName(
                "spline", spline(linearLogI(Math.log(Ii), 0.0, knots), knots, grid),
                "otherSpline", spline(linearLogI(Math.log(Ij), 0.0, knots), knots, grid),
                "forwardMigration", new RealParameter(new Double[] { m }));

        assertEquals(m * Ij / (GAMMA * Ii), f.getArrayValue(), 1e-9);
    }

    @Test
    public void testGrowingFocalDemeWithSharedMigrationVector() {
        // I_i(tau) = c * exp(-r tau) (growing forward in time), I_j constant.
        // beta_i = gamma + r; denominator = (gamma + r) * c * (1 - exp(-r T)) / r.
        RateShifts knots = rateShifts(KNOTS);
        RateShifts grid = rateShifts(fineGrid(1.0, 1000));
        double c = 200.0;
        double r = 3.0;
        double Ij = 40.0;
        double m = 1.5;
        double T = 1.0;

        // Flat 3-deme vector, focal deme i = 0, source j = 2: index j*(n-1) + i = 4.
        RealParameter migration = new RealParameter(new Double[] { 9.0, 9.0, 9.0, 9.0, m, 9.0 });

        IntroductionFraction f = new IntroductionFraction();
        f.initByName(
                "spline", spline(linearLogI(Math.log(c), -r, knots), knots, grid),
                "otherSpline", spline(linearLogI(Math.log(Ij), 0.0, knots), knots, grid),
                "forwardMigration", migration,
                "forwardMigrationIndices", new IntegerParameter(new Integer[] { 4 }));

        double expected = (m * Ij * T) / ((GAMMA + r) * c * (1 - Math.exp(-r * T)) / r);
        assertEquals(expected, f.getArrayValue(), 1e-4 * expected);
    }

    @Test
    public void testCustomWindow() {
        // Same growing focal deme as above, integrated over [a, b] instead of [0, T]:
        // numerator = m * I_j * (b - a); denominator = (gamma + r) * c * (e^{-ra} - e^{-rb}) / r.
        // The window ends deliberately fall between grid points.
        RateShifts knots = rateShifts(KNOTS);
        RateShifts grid = rateShifts(fineGrid(1.0, 1000));
        double c = 200.0;
        double r = 3.0;
        double Ij = 40.0;
        double m = 1.5;
        double a = 0.2345;
        double b = 0.7123;

        IntroductionFraction f = new IntroductionFraction();
        f.initByName(
                "spline", spline(linearLogI(Math.log(c), -r, knots), knots, grid),
                "otherSpline", spline(linearLogI(Math.log(Ij), 0.0, knots), knots, grid),
                "forwardMigration", new RealParameter(new Double[] { m }),
                "recentTime", new RealParameter(new Double[] { a }),
                "oldestTime", new RealParameter(new Double[] { b }));

        double expected = (m * Ij * (b - a))
                / ((GAMMA + r) * c * (Math.exp(-r * a) - Math.exp(-r * b)) / r);
        assertEquals(expected, f.getArrayValue(), 1e-4 * expected);
    }

    @Test
    public void testDefaultWindowEqualsKnotRange() {
        RateShifts knots = rateShifts(KNOTS);
        RateShifts grid = rateShifts(fineGrid(1.0, 100));
        Spline focal = spline(linearLogI(Math.log(200.0), -3.0, knots), knots, grid);
        Spline other = spline(linearLogI(Math.log(40.0), 0.0, knots), knots, grid);
        RealParameter migration = new RealParameter(new Double[] { 1.5 });

        IntroductionFraction byDefault = new IntroductionFraction();
        byDefault.initByName("spline", focal, "otherSpline", other, "forwardMigration", migration);
        IntroductionFraction explicit = new IntroductionFraction();
        explicit.initByName("spline", focal, "otherSpline", other, "forwardMigration", migration,
                "recentTime", new RealParameter(new Double[] { 0.0 }),
                "oldestTime", new RealParameter(new Double[] { 1.0 }));

        assertEquals(byDefault.getArrayValue(), explicit.getArrayValue(), 1e-12);
    }

    @Test
    public void testInvalidWindowIsRejected() {
        RateShifts knots = rateShifts(KNOTS);
        RateShifts grid = rateShifts(fineGrid(1.0, 100));
        Spline focal = spline(linearLogI(Math.log(100.0), 0.0, knots), knots, grid);
        Spline other = spline(linearLogI(Math.log(50.0), 0.0, knots), knots, grid);
        RealParameter migration = new RealParameter(new Double[] { 2.0 });

        // initByName wraps the IllegalArgumentException from initAndValidate in a RuntimeException.
        // recentTime must be smaller than oldestTime
        assertThrows(RuntimeException.class, () -> new IntroductionFraction().initByName(
                "spline", focal, "otherSpline", other, "forwardMigration", migration,
                "recentTime", new RealParameter(new Double[] { 0.6 }),
                "oldestTime", new RealParameter(new Double[] { 0.4 })));
        // window must lie within the knots [0, 1]
        assertThrows(RuntimeException.class, () -> new IntroductionFraction().initByName(
                "spline", focal, "otherSpline", other, "forwardMigration", migration,
                "oldestTime", new RealParameter(new Double[] { 1.5 })));
    }

    @Test
    public void testLogTracksParameterChangesOutsideModelGraph() {
        // Used as a plain <log>, no calculation-node machinery invalidates the cache,
        // so log() must reflect the current migration rate on every call.
        RateShifts knots = rateShifts(KNOTS);
        RateShifts grid = rateShifts(fineGrid(1.0, 100));
        RealParameter migration = new RealParameter(new Double[] { 2.0 });

        IntroductionFraction f = new IntroductionFraction();
        f.initByName(
                "spline", spline(linearLogI(Math.log(100.0), 0.0, knots), knots, grid),
                "otherSpline", spline(linearLogI(Math.log(50.0), 0.0, knots), knots, grid),
                "forwardMigration", migration);

        assertEquals(0.1, f.getArrayValue(), 1e-9);
        migration.setValue(0, 4.0);
        java.io.ByteArrayOutputStream bytes = new java.io.ByteArrayOutputStream();
        f.log(0, new java.io.PrintStream(bytes));
        assertEquals(0.2, Double.parseDouble(bytes.toString().trim()), 1e-9);
    }

    @Test
    public void testWorksAsPriorArgument() {
        RateShifts knots = rateShifts(KNOTS);
        RateShifts grid = rateShifts(fineGrid(1.0, 100));

        IntroductionFraction f = new IntroductionFraction();
        f.initByName(
                "spline", spline(linearLogI(Math.log(100.0), 0.0, knots), knots, grid),
                "otherSpline", spline(linearLogI(Math.log(50.0), 0.0, knots), knots, grid),
                "forwardMigration", new RealParameter(new Double[] { 2.0 }));

        Beta beta = new Beta();
        beta.initByName(
                "alpha", new RealParameter(new Double[] { 10.0 }),
                "beta", new RealParameter(new Double[] { 90.0 }));
        Prior prior = new Prior();
        prior.initByName("x", f, "distr", beta);

        // P = 0.1 for these values
        assertEquals(beta.logDensity(0.1), prior.calculateLogP(), 1e-9);
    }
}
