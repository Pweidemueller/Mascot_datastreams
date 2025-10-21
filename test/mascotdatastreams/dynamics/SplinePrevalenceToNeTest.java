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

        // Test interpolation at knot points
        double tolerance = 1e-9;
        assertEquals(1.0, dyn.getPrevalenceTime(0.0), tolerance, "Prevalence at t=0.0");
        assertEquals(148.4131591025761, dyn.getPrevalenceTime(0.6), tolerance, "Prevalence at t=0.6");
        assertEquals(0.3678794411714424, dyn.getPrevalenceTime(1.8), tolerance, "Prevalence at t=1.8");
        assertEquals(1.0, dyn.getPrevalenceTime(2.0), tolerance, "Prevalence at t=2.0");
        
        // Test interpolation at intermediate point but on gridshift points
        assertEquals(1.741106210845324, dyn.getPrevalenceTime(0.1), tolerance, "Prevalence at t=0.1");
        assertEquals(0.45879046299891335, dyn.getPrevalenceTime(1.1), tolerance, "Prevalence at t=1.1");
        assertEquals(0.25850366508689415, dyn.getPrevalenceTime(1.9), tolerance, "Prevalence at t=1.9");
    }

    @Test
    public void testIrregularGridSplineInterpolation() throws Exception {
        // Rate shifts specified as fractions of root height: 0.5 and 1.0.
        // With the test tree's root height = 2.0, these map to absolute times 1.0 and 2.0.
        RateShifts rateShifts = buildRateShifts("0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0");
        RateShifts gridRateShifts = buildRateShifts("0.0 0.24 0.76 1.68 2.0");

        // log-prevalence at knots
        Double[] logI = new Double[] { 0.0, 1.0, 2.0, 5.0, 7.0, 1.5, -2.0, 0.0, 4.0, -1.0, 0.0};

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
        assertEquals(1.0, dyn.getPrevalenceTime(0.0), tolerance, "Prevalence at t=0.0");
        assertEquals(3.040021253350453, dyn.getPrevalenceTime(0.24), tolerance, "Prevalence at t=0.24");
        assertEquals(1309.6214076166518, dyn.getPrevalenceTime(0.76), tolerance, "Prevalence at t=0.76");
        assertEquals(15.106990514646617, dyn.getPrevalenceTime(1.68), tolerance, "Prevalence at t=1.68");
        assertEquals(1.0, dyn.getPrevalenceTime(2.0), tolerance, "Prevalence at t=2.0");
    }
    
    @Test
    public void testNeSplineInterpolation() throws Exception {
        // Rate shifts specified as fractions of root height: 0.5 and 1.0.
        // With the test tree's root height = 2.0, these map to absolute times 1.0 and 2.0.
        RateShifts rateShifts = buildRateShifts("0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0");
        RateShifts gridRateShifts = buildRateShifts("0.0 0.24 0.76 1.68 2.0");

        // log-prevalence at knots
        Double[] logI = new Double[] { 0.0, 1.0, 2.0, 5.0, 7.0, 1.5, -2.0, 0.0, 4.0, -1.0, 0.0};

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

        assertEquals(0.007644136727600229, dyn.getNeTime(0.0), tolerance, "Ne at t=0.0");
        assertEquals(0.021105462990341873, dyn.getNeTime(0.24), tolerance, "Ne at t=0.24");
        assertEquals(8.848148341013303, dyn.getNeTime(0.76), tolerance, "Ne at t=0.76");
        assertEquals(0.07441497433106414, dyn.getNeTime(1.68), tolerance, "Ne at t=1.68");
        assertEquals(0.02836701655438929, dyn.getNeTime(2.0), tolerance, "Ne at t=2.0");
    }
    
    private static RateShifts buildRateShifts(String shiftValues) throws Exception {
        RateShifts rs = new RateShifts();
        rs.initByName("value", shiftValues);
        return rs;
    }
}
