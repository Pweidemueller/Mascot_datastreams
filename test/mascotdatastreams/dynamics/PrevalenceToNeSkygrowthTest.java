package mascotdatastreams.dynamics;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import mascot.parameterdynamics.Skygrowth;

/**
 * Tests for PrevalenceToNeSkygrowth.getNeTime().
 * We construct a simple set of rate shifts and log-prevalence values, then
 * verify Ne(t) matches manual calculations, including behavior after the last
 * shift and with custom coalescent scaling.
 */
public class PrevalenceToNeSkygrowthTest {

    @Test
    public void testNeAcrossIntervals_defaultC() throws Exception {
        // Rate shifts specified as fractions of root height: 0.5 and 1.0.
        // With the test tree's root height = 2.0, these map to absolute times 1.0 and 2.0.
        Object rateShifts = buildRateShifts("0.5 1.0");

        // log-prevalence at control points: [ln 10, ln 20, ln 5]
        Double[] logI = new Double[] { Math.log(10.0), Math.log(20.0), Math.log(5.0) };

        // gamma (uninfectious rate) = 1.0, coalescentScale uses default c=2.0
        RealParameter gamma = new RealParameter(new Double[] { 1.0 });

        // Build prevalence dynamics (Mascot Skygrowth) and wire into PrevalenceToNeSkygrowth
        Skygrowth prev = new Skygrowth();
        prev.initByName(
                "logNe", new RealParameter(logI),
                "rateShifts", rateShifts
        );
        prev.initAndValidate();

        PrevalenceToNeSkygrowth dyn = new PrevalenceToNeSkygrowth();
        dyn.initByName(
                "prevalence", prev,
                "uninfectiousRate", gamma
        );
        dyn.initAndValidate();

        // Pre-computed growth (forward-time) from implementation:
        // growth[0] = ln(10) - ln(20) over dt=1.0 => ln(0.5) ~ -0.69314718056
        // growth[1] = ln(20) - ln(5)  over dt=1.0 => ln(4)   ~  1.38629436112
        // For t=0.5 (interval 0):
        //   logI(t) = ln(10) - (-0.6931)*0.5 = ln(10) + 0.3466 -> I(t) ~ 14.14213562
        //   transmission = g_fwd + gamma = -0.6931 + 1.0 = 0.30685281988
        //   Ne = I / (2 * transmission)
        double t1 = 0.5;
        double ne1 = dyn.getNeTime(t1);

        double growth1 = (Math.log(10.0)-Math.log(20.0))/1.0;
        double I1 = 10.0 / Math.exp(growth1*t1); // ~14.14213562
        double transmission1 = growth1 + gamma.getArrayValue(); // ~0.30685281944
        double expected1 = I1 / (2.0 * transmission1);
        assertEquals(expected1, ne1, 1e-9, "Ne(t) at t=0.5 incorrect");

        // For t=1.5 (interval 1):
        //   logI(t) = ln(20) - 1.3863*0.5 => I(t)=10
        //   transmission = 1.3863 + 1.0 = 2.38629436112
        //   Ne = 10 / (2 * 2.38629...)
        double t2 = 1.5;
        double ne2 = dyn.getNeTime(t2);
        double growth2 = (Math.log(20.0)-Math.log(5.0))/1.0;
        double I2 = 20.0 / Math.exp(growth2*0.5); 
        double transmission2 = growth2 + gamma.getArrayValue();
        double expected2 = I2 / (2.0 * transmission2);
        assertEquals(expected2, ne2, 1e-9, "Ne(t) at t=1.5 incorrect");

        // For t=3.0 (after last shift):
        //   logI(t) = ln(5) (last control point)
        //   g_fwd is 0
        //   transmission = 0 + 1.0; I=5; Ne = 5 / (2 * transmission)
        double t3 = 3.0;
        double ne3 = dyn.getNeTime(t3);
        double growth3 = 0.0;
        double I3 = 5.0 ; 
        double transmission3 = growth3 + gamma.getArrayValue();
        double expected3 = I3 / (2.0 * transmission3);
        assertEquals(expected3, ne3, 1e-9, "Ne(t) after last shift incorrect");
        
        // For t=0.0 (at the most recent sample):
        //   logI(t) = ln(10) (first control point)
        //   growth[0] = ln(10) - ln(20) over dt=1.0 => ln(0.5) ~ -0.69314718056
        //   transmission = -0.6931 + 1.0; I=10; Ne = 10 / (2 * transmission)
        double t4 = 0.0;
        double ne4 = dyn.getNeTime(t4);
        double growth4 = (Math.log(10.0)-Math.log(20.0))/1.0;
        double I4 = 10.0 ; 
        double transmission4 = growth4 + gamma.getArrayValue();
        double expected4 = I4 / (2.0 * transmission4);
        assertEquals(expected4, ne4, 1e-9, "Ne(t) at t=0.0 incorrect");
    }

    @Test
    public void testNeWithCustomCoalescentScale() throws Exception {
        // Fractions of root height (2.0): 0.5->1.0, 1.0->2.0
        Object rateShifts = buildRateShifts("0.5 1.0");
        Double[] logI = new Double[] { Math.log(10.0), Math.log(20.0), Math.log(5.0) };
        RealParameter gamma = new RealParameter(new Double[] { 1.0 });
        RealParameter c = new RealParameter(new Double[] { 4.0 }); // double the denominator vs default

        Skygrowth prev = new Skygrowth();
        prev.initByName(
                "logNe", new RealParameter(logI),
                "rateShifts", rateShifts
        );
        prev.initAndValidate();

        PrevalenceToNeSkygrowth dyn = new PrevalenceToNeSkygrowth();
        dyn.initByName(
                "prevalence", prev,
                "uninfectiousRate", gamma,
                "coalescentScale", c
        );
        dyn.initAndValidate();

        // Compare to default-c case at t=1.5, where I=10 and g_fwd=1.38629
        // transmission = 1.38629 + 1.0 = 2.38629
        // default Ne_default = 10 / (2 * 2.38629)
        // with c=4.0, Ne_c4 = 10 / (4 * 2.38629) = 0.5 * Ne_default
        double t = 1.5;
        double ne_c4 = dyn.getNeTime(t);

        // Compute expected using same numbers
        double growth = (Math.log(20.0)-Math.log(5.0))/1.0;
        double I = 20.0 / Math.exp(growth*0.5); 
        double transmission = growth + gamma.getArrayValue();
        double expected_c4 = I / (c.getArrayValue() * transmission);
        assertEquals(expected_c4, ne_c4, 1e-9, "Ne(t) with custom c incorrect");
    }

    // Build a minimal tree (root height = 2.0) and RateShifts instance initialized with fractional breakpoints
    private static Object buildRateShifts(String shiftValues) throws Exception {
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

        Class<?> clsRateShifts = Class.forName("mascot.dynamics.RateShifts");
        Object rs = clsRateShifts.getDeclaredConstructor().newInstance();
        clsRateShifts.getMethod("initByName", Object[].class)
                .invoke(rs, new Object[] { new Object[] { "tree", tree, "value", shiftValues } });
        return rs;
    }
}
