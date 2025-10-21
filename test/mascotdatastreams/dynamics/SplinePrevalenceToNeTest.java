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
        RateShifts rateShifts = buildRateShifts("0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0");
        RateShifts gridRateShifts = buildRateShifts("0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 0.99 1.0 1.05 1.1");

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
        RateShifts rateShifts = buildRateShifts("0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0");
        RateShifts gridRateShifts = buildRateShifts("0.0 0.12 0.38 0.84 1.0");

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
        RateShifts rateShifts = buildRateShifts("0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0");
        RateShifts gridRateShifts = buildRateShifts("0.0 0.12 0.38 0.84 1.0");

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
    
    // Build a minimal tree (root height = 2.0) and RateShifts instance initialized with fractional breakpoints
    private static RateShifts buildRateShifts(String shiftValues) throws Exception {
        // Minimal BEAST setup copied from existing tests to satisfy RateShifts
        Sequence s1 = new Sequence();
        Sequence s2 = new Sequence();
        Sequence s3 = new Sequence();
        Sequence s4 = new Sequence();
        s1.initByName("taxon", "a1", "value", "?");
        s2.initByName("taxon", "b1", "value", "?");
        s3.initByName("taxon", "a2", "value", "?");
        s4.initByName("taxon", "b2", "value", "?");
        Alignment alignment = new Alignment();
        alignment.initByName("sequence", s1, "sequence", s2, "sequence", s3, "sequence", s4);

        TaxonSet taxa = new TaxonSet();
        taxa.initByName("alignment", alignment);

        TraitSet traitSet = new TraitSet();
        traitSet.initByName("value", "a1=Deme1,a2=Deme1,b1=Deme2,b2=Deme2", "traitname", "type", "taxa", taxa);

        // 4-tip ultrametric tree; branch lengths sum to 2.0 from root to tips (root height = 2.0)
        // Left clade: (a1:1.0,a2:1.0):1.0 -> total to tips = 2.0
        // Right clade: (b1:1.0,b2:1.0):1.0 -> total to tips = 2.0
        Tree tree = new TreeParser("((a1:1.0,a2:1.0):1.0,(b1:1.0,b2:1.0):1.0);");
        tree.initByName("taxonset", taxa, "trait", traitSet);

        RateShifts rs = new RateShifts();
        rs.initByName("tree", tree, "value", shiftValues);
        return rs;
    }
}
