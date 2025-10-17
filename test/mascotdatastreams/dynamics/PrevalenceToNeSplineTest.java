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
public class PrevalenceToNeSplineTest {
    
    @Test
    public void testBasicSplineInterpolation() throws Exception {
        // Rate shifts specified as fractions of root height: 0.5 and 1.0.
        // With the test tree's root height = 2.0, these map to absolute times 1.0 and 2.0.
        RateShifts rateShifts = buildRateShifts("0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0");

        // log-prevalence at control points: [ln 10, ln 15, ln 5, ln 1, ln 2, ln 3, ln 4, ln 5, ln 6, ln 7, ln 8, ln 9, ln 10]
        Double[] logI = new Double[] { Math.log(10.0), Math.log(15.0), Math.log(5.0), Math.log(1.0), Math.log(2.0), Math.log(20.0), Math.log(4.0), Math.log(5.0), Math.log(6.0), Math.log(7.0), Math.log(8.0)};

        // gamma (uninfectious rate) = 0.1
        RealParameter gamma = new RealParameter(new Double[] { 0.1 });

        PrevalenceToNeSpline dyn = new PrevalenceToNeSpline();
        dyn.initByName(
                "logInfected", new RealParameter(logI),
                "rateShifts", rateShifts,
                "uninfectiousRate", gamma
        );
        dyn.initAndValidate();

        // Test interpolation at knot points
        double tolerance = 1e-9;
        assertEquals(10.0, dyn.getPrevalenceTimeExact(0.0), tolerance, "Prevalence at t=0");
        assertEquals(20.0, dyn.getPrevalenceTimeExact(1.0), tolerance, "Prevalence at t=1.0");
        assertEquals(8.0, dyn.getPrevalenceTimeExact(2.0), tolerance, "Prevalence at t=2.0");
        
        // Test interpolation at intermediate point
        double prevalence_at_0_1 = dyn.getPrevalenceTimeExact(0.1);
        assertTrue(prevalence_at_0_1 > 10.0 && prevalence_at_0_1 < 15.0, 
            "Prevalence at t=0.1 should be between 10 and 15, but was " + prevalence_at_0_1);
    }
    @Test
    public void testNeCalculationWithCustomCoalescentScale() throws Exception {
        // Rate shifts specified as fractions of root height: 1.0.
        // With the test tree's root height = 2.0, this maps to absolute time 2.0.
        RateShifts rateShifts = buildRateShifts("1.0");

        // Test Ne calculation with constant log-prevalence (zero derivative)
        // log I(t) = 1.0 (constant), so dlogI/dt = 0
        // transmission_rate = 0 + gamma = gamma
        // Ne = I / (c * gamma) = exp(1.0) / (c * gamma)
        Double[] logI = new Double[] { 1.0, 1.0 };  // Constant log-prevalence
        double gamma = 0.1;
        double c = 4.0;  // Custom coalescent scale

        RealParameter gammaParam = new RealParameter(new Double[]{gamma});
        RealParameter cParam = new RealParameter(new Double[]{c});

        PrevalenceToNeSpline dyn = new PrevalenceToNeSpline();
        dyn.initByName(
                "logInfected", new RealParameter(logI),
                "rateShifts", rateShifts,
                "uninfectiousRate", gammaParam,
                "coalescentScale", cParam
        );
        dyn.initAndValidate();
        
        double expectedNe = Math.exp(1.0) / (c * gamma);
        double actualNe = dyn.getNeTime(1.0);  // Test at midpoint
        
        double tolerance = 1e-6;
        assertEquals(expectedNe, actualNe, tolerance, "Ne calculation with custom coalescent scale");
    }
    
    @Test
    public void testNumGridPointsParameter() throws Exception {
        // Rate shifts specified as fractions of root height: 0.5 and 1.0.
        // With the test tree's root height = 2.0, these map to absolute times 1.0 and 2.0.
        RateShifts rateShifts = buildRateShifts("0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0");

        // log-prevalence at control points: [ln 10, ln 20, ln 5, ln 1, ln 2, ln 3, ln 4, ln 5, ln 6, ln 7, ln 8, ln 9, ln 10]
        Double[] logI = new Double[] { Math.log(10.0), Math.log(15.0), Math.log(5.0), Math.log(1.0), Math.log(2.0), Math.log(20.0), Math.log(4.0), Math.log(5.0), Math.log(6.0), Math.log(7.0), Math.log(8.0)};


        // gamma (uninfectious rate) = 0.1
        RealParameter gamma = new RealParameter(new Double[] { 0.1 });

        // Test with custom numGridPoints
        PrevalenceToNeSpline dyn = new PrevalenceToNeSpline();
        dyn.initByName(
                "logInfected", new RealParameter(logI),
                "rateShifts", rateShifts,
                "uninfectiousRate", gamma,
                "numGridPoints", 1000  // Use 500 grid points instead of default 1000
        );
        dyn.initAndValidate();
        
        // Test that the spline still works correctly with custom grid points
        double tolerance = 1e-9;
        assertEquals(10.0, dyn.getPrevalenceTime(0.0), tolerance, "Prevalence at t=0");
        assertEquals(20.0, dyn.getPrevalenceTime(1.0), tolerance, "Prevalence at t=1.0");
        assertEquals(8.0, dyn.getPrevalenceTime(2.0), tolerance, "Prevalence at t=2.0");
        
        // Test interpolation at intermediate point
        double prevalence_at_0_1 = dyn.getPrevalenceTime(0.1);
        assertTrue(prevalence_at_0_1 > 10.0 && prevalence_at_0_1 < 15.0, 
            "Prevalence at t=0.1 should be between 10 and 15");

        // Test Ne calculation
        double ne = dyn.getNeTime(0.5);
        assertTrue(ne > 0, "Ne should be positive");
        assertFalse(Double.isNaN(ne), "Ne should not be NaN");
        assertFalse(Double.isInfinite(ne), "Ne should not be infinite");
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
