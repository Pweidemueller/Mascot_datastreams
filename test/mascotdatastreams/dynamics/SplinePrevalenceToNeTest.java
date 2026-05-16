package mascotdatastreams.dynamics;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.assertFalse;

import java.util.Arrays;

import org.junit.jupiter.api.Test;

import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import mascot.dynamics.RateShifts;

/**
 * Unit tests for PrevalenceToNeSpline class.
 * 
 * Tests spline interpolation accuracy, Ne calculation with known prevalence trajectories,
 * boundary conditions, and parameter validation. Uses proper tree setup and rate shifts
 * following the pattern from PrevalenceToNeSkygrowthTest.
 */
public class SplinePrevalenceToNeTest {
    
    @Test
    public void testBasicSplineInterpolation() throws Exception {
        // Rate shifts specified as fractions of root height: 0.5 and 1.0.
        // With the test tree's root height = 2.0, these map to absolute times 1.0 and 2.0.
    	RateShifts rateShifts = buildRateShifts("0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0");
        RateShifts gridRateShifts = buildRateShifts("0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0");

        Double[] logI = new Double[] { 0.0, 1.0, 2.0, 5.0, 10.0, 2.5, -2.0, 0.0, 1.0, -1.0, 0.0};

        // gamma (uninfectious rate)
        RealParameter gamma = new RealParameter(new Double[] { 75.0 });

        // Create the Spline first
        Spline spline = new Spline();
        spline.initByName(
                "logInfected", new RealParameter(logI),
                "rateShifts", rateShifts,
                "gridRateShifts", gridRateShifts,
                "uninfectiousRate", gamma
        );
        spline.initAndValidate();

        // Create SplinePrevalenceToNe with the spline
        SplinePrevalenceToNe dyn = new SplinePrevalenceToNe();
        dyn.initByName(
                "spline", spline,
                "coalescentScale", new RealParameter(new Double[] { 2.0 })
        );
        dyn.initAndValidate();

        // Test interpolation at knot points
        double tolerance = 1e-9;
        assertEquals(0.0, Math.log(dyn.getPrevalenceTime(0.0)), tolerance, "Prevalence at t=0.0");
        assertEquals(5.0, Math.log(dyn.getPrevalenceTime(0.6)), tolerance, "Prevalence at t=0.6");
        assertEquals(10.0, Math.log(dyn.getPrevalenceTime(0.8)), tolerance, "Prevalence at t=0.8");
        assertEquals(2.5, Math.log(dyn.getPrevalenceTime(1.0)), tolerance, "Prevalence at t=1.0");
        assertEquals(-1.0, Math.log(dyn.getPrevalenceTime(1.8)), tolerance, "Prevalence at t=1.8");
        assertEquals(0.0, Math.log(dyn.getPrevalenceTime(2.0)), tolerance, "Prevalence at t=2.0");
        
        // Test interpolation at intermediate point but on gridshift points
        assertEquals(1.6729764928095872, dyn.getPrevalenceTime(0.1), tolerance, "Prevalence at t=0.1");
        assertEquals(2.654981149109626, dyn.getPrevalenceTime(1.5), tolerance, "Prevalence at t=1.5");
        
        // Test interpolation at times not on grid
        assertEquals(1.0528075333474662, dyn.getPrevalenceTime(0.01), tolerance, "Prevalence at t=0.01");
        assertEquals(1.2934359252818004, dyn.getPrevalenceTime(0.05), tolerance, "Prevalence at t=0.05");
        assertEquals(1.5890620458330644, dyn.getPrevalenceTime(0.09), tolerance, "Prevalence at t=0.09");
        assertEquals(0.4587726242235013, dyn.getPrevalenceTime(1.91), tolerance, "Prevalence at t=1.91");
        assertEquals(0.6486323549478357, dyn.getPrevalenceTime(1.95), tolerance, "Prevalence at t=1.95");
        assertEquals(0.8410066760378204, dyn.getPrevalenceTime(1.98), tolerance, "Prevalence at t=1.98");
        
        // Test at times outside the knot intervals
        assertEquals(1.0, dyn.getPrevalenceTime(0.0), tolerance, "Prevalence at t=-0.1");
        assertEquals(1.0, dyn.getPrevalenceTime(2.05), tolerance, "Prevalence at t=2.05");
        
    }

   @Test
   public void testIrregularGridSplineInterpolation() throws Exception {
        // Rate shifts specified as fractions of root height: 0.5 and 1.0.
        // With the test tree's root height = 2.0, these map to absolute times 1.0 and 2.0.
        RateShifts rateShifts = buildRateShifts("0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0");
        RateShifts gridRateShifts = buildRateShifts("0.0 0.24 0.76 1.68 2.0");

        // log-prevalence at knots
        Double[] logI = new Double[] { 0.0, 1.0, 2.0, 5.0, 10.0, 2.5, -2.0, 0.0, 1.0, -1.0, 0.0};

        // gamma (uninfectious rate)
        RealParameter gamma = new RealParameter(new Double[] { 75.0 });

        // Create the Spline first
        Spline spline = new Spline();
        spline.initByName(
                "logInfected", new RealParameter(logI),
                "rateShifts", rateShifts,
                "gridRateShifts", gridRateShifts,
                "uninfectiousRate", gamma
        );
        spline.initAndValidate();

        // Create SplinePrevalenceToNe with the spline
        SplinePrevalenceToNe dyn = new SplinePrevalenceToNe();
        dyn.initByName(
                "spline", spline,
                "coalescentScale", new RealParameter(new Double[] { 2.0 })
        );
        dyn.initAndValidate();
        
        // Test interpolation at intermediate point but on gridshift points
        double tolerance = 1e-9;
        assertEquals(3.2586285123931096, dyn.getPrevalenceTime(0.24), tolerance, "Prevalence at t=0.24");
        assertEquals(18546.878553616978, dyn.getPrevalenceTime(0.76), tolerance, "Prevalence at t=0.76");
        assertEquals(1.2441976859534918, dyn.getPrevalenceTime(1.68), tolerance, "Prevalence at t=1.68");
        
        // Test interpolation at times not on grid
        assertEquals(1.0504525803565143, dyn.getPrevalenceTime(0.01), tolerance, "Prevalence at t=0.01");
        assertEquals(1.2790345059335768, dyn.getPrevalenceTime(0.05), tolerance, "Prevalence at t=0.05");
        assertEquals(1.557356607961807,  dyn.getPrevalenceTime(0.09), tolerance, "Prevalence at t=0.09");
        assertEquals(2.0924099822312656, dyn.getPrevalenceTime(0.15), tolerance, "Prevalence at t=0.15");
        assertEquals(208.1789262746285,  dyn.getPrevalenceTime(0.49), tolerance, "Prevalence at t=0.49");
        assertEquals(5.366400513902776,  dyn.getPrevalenceTime(0.27), tolerance, "Prevalence at t=0.27");
   }
   
   @Test
   public void testNeSplineInterpolation() throws Exception {
        // Rate shifts specified as fractions of root height: 0.5 and 1.0.
        // With the test tree's root height = 2.0, these map to absolute times 1.0 and 2.0.
        RateShifts rateShifts = buildRateShifts("0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0");
        RateShifts gridRateShifts = buildRateShifts("0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0");

        // log-prevalence at knots
        Double[] logI = new Double[] { 0.0, 1.0, 2.0, 5.0, 10.0, 2.5, -2.0, 0.0, 1.0, -1.0, 0.0};

        // gamma (uninfectious rate)
        RealParameter gamma = new RealParameter(new Double[] { 75.0 });

        // Create the Spline first
        Spline spline = new Spline();
        spline.initByName(
                "logInfected", new RealParameter(logI),
                "rateShifts", rateShifts,
                "gridRateShifts", gridRateShifts,
                "uninfectiousRate", gamma
        );
        spline.initAndValidate();

        // Create SplinePrevalenceToNe with the spline
        SplinePrevalenceToNe dyn = new SplinePrevalenceToNe();
        dyn.initByName(
                "spline", spline,
                "coalescentScale", new RealParameter(new Double[] { 2.0 })
        );
        dyn.initAndValidate();
        
        // Test interpolation at intermediate point but on gridshift points
        double tolerance = 1e-9;

        assertEquals(0.007162782463310597, dyn.getNeTime(0.0), tolerance, "Ne at t=0.0");
        assertEquals(0.011958148340766807, dyn.getNeTime(0.1), tolerance, "Ne at t=0.1");
        assertEquals(0.01914630507153202, dyn.getNeTime(1.5), tolerance, "Ne at t=1.5");
        assertEquals(0.007677784323246716, dyn.getNeTime(2.0), tolerance, "Ne at t=2.0");
        
        // outside grid points
        assertEquals(0.12153367392149209,  dyn.getNeTime(0.48),  tolerance, "Ne at t=0.48");          
        assertEquals(0.0024361235349888563, dyn.getNeTime(1.82),  tolerance, "Ne at t=1.82");
        
        // outside the knot intervals
        assertEquals(0.007162782463310597, dyn.getNeTime(-0.1),  tolerance, "Ne at t=-0.1");          
        assertEquals(0.007677784323246716, dyn.getNeTime(2.05),  tolerance, "Ne at t=2.05");     

        // Apply NeScaler and ensure Ne(t) values scale multiplicatively.
        SplinePrevalenceToNe dynScaled = new SplinePrevalenceToNe();
        dynScaled.initByName(
                "spline", spline,
                "coalescentScale", new RealParameter(new Double[] { 2.0 }),
                "NeScaler", new RealParameter(new Double[] { 2.0 })
        );
        dynScaled.initAndValidate();

        assertEquals(0.007162782463310597 * 2.0, dynScaled.getNeTime(0.0), tolerance, "NeScaler at t=0.0");
        assertEquals(0.011958148340766807 * 2.0, dynScaled.getNeTime(0.1), tolerance, "NeScaler at t=0.1");
        assertEquals(0.01914630507153202 * 2.0, dynScaled.getNeTime(1.5), tolerance, "NeScaler at t=1.5");
        assertEquals(0.007677784323246716 * 2.0, dynScaled.getNeTime(2.0), tolerance, "NeScaler at t=2.0");
   }
  
    @Test
    public void testLocalTransmissionRateAndNe() throws Exception {
        // Two-deme setup. Deme i is "self", deme j is "other". A single forward
        // migration rate m_{j -> i} subtracts from i's aggregate transmission rate.
        RateShifts rateShifts = buildRateShifts("0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0");
        RateShifts gridRateShifts = buildRateShifts("0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0");

        Double[] logI_i = new Double[] { 0.0, 1.0, 2.0, 5.0, 10.0, 2.5, -2.0, 0.0, 1.0, -1.0, 0.0 };
        Double[] logI_j = new Double[] { 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.3, 1.1, 0.9, 0.7, 0.5 };

        RealParameter gamma = new RealParameter(new Double[] { 75.0 });

        Spline splineI = new Spline();
        splineI.initByName(
                "logInfected", new RealParameter(logI_i),
                "rateShifts", rateShifts,
                "gridRateShifts", gridRateShifts,
                "uninfectiousRate", gamma);
        splineI.initAndValidate();

        Spline splineJ = new Spline();
        splineJ.initByName(
                "logInfected", new RealParameter(logI_j),
                "rateShifts", rateShifts,
                "gridRateShifts", gridRateShifts,
                "uninfectiousRate", gamma);
        splineJ.initAndValidate();

        // Baseline (no migration subtraction) so we can compute the expected
        // aggregate transmission rate via the same code path as the production class.
        SplinePrevalenceToNe baseline = new SplinePrevalenceToNe();
        baseline.initByName(
                "spline", splineI,
                "coalescentScale", new RealParameter(new Double[] { 2.0 }));
        baseline.initAndValidate();

        double m_ji = 0.3;
        SplinePrevalenceToNe withMig = new SplinePrevalenceToNe();
        withMig.initByName(
                "spline", splineI,
                "coalescentScale", new RealParameter(new Double[] { 2.0 }),
                "otherSpline", Arrays.asList(splineJ),
                "incomingForwardMigration", new RealParameter(new Double[] { m_ji }));
        withMig.initAndValidate();

        double tolerance = 1e-9;
        double[] testTimes = { 0.0, 0.1, 0.55, 1.0, 1.5, 2.0 };
        for (double t : testTimes) {
            double I_i = splineI.getPrevalence(t);
            double I_j = splineJ.getPrevalence(t);
            double betaTotal = splineI.getTransmissionRate(t);
            double expectedMig = m_ji * I_j / I_i;
            double expectedLocal = betaTotal - expectedMig;

            assertEquals(expectedMig, withMig.getMigrationTransmissionRate(t), tolerance,
                    "migration rate at t=" + t);
            assertEquals(expectedLocal, withMig.getLocalTransmissionRate(t), tolerance,
                    "local rate at t=" + t);

            // Ne should use the local rate when migration info is wired.
            double expectedNe = I_i / (2.0 * expectedLocal);
            assertEquals(expectedNe, withMig.getNeTime(t), tolerance, "Ne at t=" + t);

            // Baseline still uses the aggregate rate.
            assertEquals(betaTotal, baseline.getLocalTransmissionRate(t), tolerance,
                    "baseline local == total at t=" + t);
            assertEquals(0.0, baseline.getMigrationTransmissionRate(t), tolerance,
                    "baseline migration is zero at t=" + t);
            assertEquals(I_i / (2.0 * betaTotal), baseline.getNeTime(t), tolerance,
                    "baseline Ne uses total beta at t=" + t);
        }
    }

    private static RateShifts buildRateShifts(String shiftValues) throws Exception {
        RateShifts rs = new RateShifts();
        rs.initByName("value", shiftValues);
        return rs;
    }
}
