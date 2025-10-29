package mascotdatastreams.dynamics;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.assertFalse;

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

        // Create the NotAKnotSpline first
        NotAKnotSpline spline = new NotAKnotSpline();
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
        assertEquals(dyn.getPrevalenceTime(0.0), dyn.getPrevalenceTime(0.01), tolerance, "Prevalence at t=0.01");
        assertEquals(dyn.getPrevalenceTime(0.0), dyn.getPrevalenceTime(0.05), tolerance, "Prevalence at t=0.05");
        assertEquals(dyn.getPrevalenceTime(0.1), dyn.getPrevalenceTime(0.09), tolerance, "Prevalence at t=0.09");
        
        assertEquals(dyn.getPrevalenceTime(1.9), dyn.getPrevalenceTime(1.91), tolerance, "Prevalence at t=1.9");
        assertEquals(dyn.getPrevalenceTime(1.9), dyn.getPrevalenceTime(1.95), tolerance, "Prevalence at t=1.95");
        assertEquals(dyn.getPrevalenceTime(2.0), dyn.getPrevalenceTime(1.98), tolerance, "Prevalence at t=1.98");
        
        // Test at times outside the knot intervals
        assertEquals(dyn.getPrevalenceTime(-0.1), dyn.getPrevalenceTime(0.0), tolerance, "Prevalence at t=-0.1");
        assertEquals(dyn.getPrevalenceTime(2.0), dyn.getPrevalenceTime(2.05), tolerance, "Prevalence at t=2.05");
        
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

        // Create the NotAKnotSpline first
        NotAKnotSpline spline = new NotAKnotSpline();
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
        assertEquals(dyn.getPrevalenceTime(0.0), dyn.getPrevalenceTime(0.01), tolerance, "Prevalence at t=0.01");
        assertEquals(dyn.getPrevalenceTime(0.0), dyn.getPrevalenceTime(0.05), tolerance, "Prevalence at t=0.05");
        assertEquals(dyn.getPrevalenceTime(0.1), dyn.getPrevalenceTime(0.09), tolerance, "Prevalence at t=0.09");
        
        assertEquals(dyn.getPrevalenceTime(0.24), dyn.getPrevalenceTime(0.15), tolerance, "Prevalence at t=1.9");
        assertEquals(dyn.getPrevalenceTime(0.24), dyn.getPrevalenceTime(0.49), tolerance, "Prevalence at t=1.95");
        assertEquals(dyn.getPrevalenceTime(0.24), dyn.getPrevalenceTime(0.27), tolerance, "Prevalence at t=1.98");
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

        // Create the NotAKnotSpline first
        NotAKnotSpline spline = new NotAKnotSpline();
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
        assertEquals(dyn.getNeTime(0.5), dyn.getNeTime(0.48), tolerance, "Ne at t=0.48");
        assertEquals(dyn.getNeTime(1.8), dyn.getNeTime(1.82), tolerance, "Ne at t=1.82");
        
        // outside the knot intervals
        assertEquals(dyn.getNeTime(-0.1), dyn.getNeTime(0.0), tolerance, "Ne at t=-0.1");
        assertEquals(dyn.getNeTime(2.0), dyn.getNeTime(2.05), tolerance, "Ne at t=2.05");
   }
  
    private static RateShifts buildRateShifts(String shiftValues) throws Exception {
        RateShifts rs = new RateShifts();
        rs.initByName("value", shiftValues);
        return rs;
    }
}
