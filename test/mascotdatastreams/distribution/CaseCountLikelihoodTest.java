package mascotdatastreams.distribution;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.condition.EnabledIf;

import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;

// The following imports target native MASCOT classes referenced in BEAST XMLs.
// They must be available on the test classpath via the MASCOT plugin/library.
// If they are not present, this test can be skipped using @EnabledIf below.
// import mascot.util.InitializedNeDynamicsList;
// import mascot.parameterdynamics.Skygrowth;
// import mascot.dynamics.RateShifts;

public class CaseCountLikelihoodTest {

    /**
     * This test uses native BEAST2/MASCOT components to build Ne dynamics:
     * - InitializedNeDynamicsList
     * - Skygrowth per deme with RateShifts and logNe skyline values
     * If these classes are not present on the classpath, the test is skipped.
     */
    @Test
    @EnabledIf("mascotClassesPresent")
    public void testTwoDemesFiveTimePointsWithSkyline() throws Exception {
        // 1) Prepare case count observations: 2 demes x 5 time points
        // Times in years (relative to most recent sample)
        Double[] times = new Double[] {0.0, 0.1, 0.11, 0.24, 0.3,
                                       0.0, 0.15, 0.16, 0.25, 0.37};
        // Trait indices: first 5 are deme 0, next 5 are deme 1
        Double[] traits = new Double[] {0.0, 0.0, 0.0, 0.0, 0.0,
                                        1.0, 1.0, 1.0, 1.0, 1.0};
        // Case counts (integers represented as doubles)
        Double[] counts = new Double[] {5.0, 7.0, 6.0, 8.0, 10.0,
                                        3.0, 4.0, 5.0, 6.0, 8.0};

        CaseCountData caseData = new CaseCountData();
        caseData.initByName(
                "caseCounts", new RealParameter(counts),
                "observationTimes", new RealParameter(times),
                "traitIndices", new RealParameter(traits)
        );

        // 2) Prepare Ne skyline with native MASCOT classes using BEAST-style initByName
        // Build a tiny tree and trait set to satisfy RateShifts' expected inputs
        Sequence s1 = new Sequence();
        Sequence s2 = new Sequence();
        Sequence s3 = new Sequence();
        Sequence s4 = new Sequence();
        s1.initByName("taxon", "a1", "value", "???");
        s2.initByName("taxon", "b1", "value", "???");
        s3.initByName("taxon", "a2", "value", "???");
        s4.initByName("taxon", "b2", "value", "???");
        Alignment alignment = new Alignment();
        alignment.initByName("sequence", s1, "sequence", s2, "sequence", s3, "sequence", s4);

        TaxonSet taxa = new TaxonSet();
        taxa.initByName("alignment", alignment);

        TraitSet traitSet = new TraitSet();
        traitSet.initByName("value", "a1=Deme1,a2=Deme1,b1=Deme2,b2=Deme2", "traitname", "type", "taxa", taxa);
        
        // 4-tip ultrametric tree with root height ~0.4
        // Left clade: (a1:0.1,a2:0.1):0.3 -> total to tips = 0.4
        // Right clade: (b1:0.2,b2:0.2):0.2 -> total to tips = 0.4
        Tree tree = new TreeParser("((a1:0.1,a2:0.1):0.3,(b1:0.2,b2:0.2):0.2);");
        tree.initByName("taxonset", taxa, "trait", traitSet);

        // Instantiate via reflection to avoid hard compile dependency if the plugin isn't present yet
        Class<?> clsRateShifts = Class.forName("mascot.dynamics.RateShifts");
        Object rateShifts1 = clsRateShifts.getDeclaredConstructor().newInstance();
        Object rateShifts2 = clsRateShifts.getDeclaredConstructor().newInstance();
        clsRateShifts.getMethod("initByName", Object[].class)
                .invoke(rateShifts1, new Object[]{new Object[]{"tree", tree, "value", "0.25 0.5 0.75 1.00"}});
        clsRateShifts.getMethod("initByName", Object[].class)
                .invoke(rateShifts2, new Object[]{new Object[]{"tree", tree, "value", "0.25 0.5 0.75 1.00"}});

        Class<?> clsSkygrowth = Class.forName("mascot.parameterdynamics.Skygrowth");
        Object skygrowth1 = clsSkygrowth.getDeclaredConstructor().newInstance();
        Object skygrowth2 = clsSkygrowth.getDeclaredConstructor().newInstance();
        // Use 3 logNe values to match 3 intervals
        RealParameter logNeDeme1 = new RealParameter(new Double[]{Math.log(5.0), Math.log(7.5), Math.log(9.0), Math.log(3.0), Math.log(1.0)});
        RealParameter logNeDeme2 = new RealParameter(new Double[]{Math.log(3.0), Math.log(4.5), Math.log(6.0), Math.log(1.0), Math.log(0.5)});
        clsSkygrowth.getMethod("initByName", Object[].class)
                .invoke(skygrowth1, new Object[]{new Object[]{"logNe", logNeDeme1, "rateShifts", rateShifts1}});
        clsSkygrowth.getMethod("initByName", Object[].class)
                .invoke(skygrowth2, new Object[]{new Object[]{"logNe", logNeDeme2, "rateShifts", rateShifts2}});

        Class<?> clsNeDynList = Class.forName("mascot.util.InitializedNeDynamicsList");
        Object neList = clsNeDynList.getDeclaredConstructor().newInstance();
        // Add both Skygrowth dynamics as children
        clsNeDynList.getMethod("initByName", Object[].class)
                .invoke(neList, new Object[]{new Object[]{"neDynamics", skygrowth1, "neDynamics", skygrowth2}});

        // 3) Prepare Gamma-Poisson distribution with overdispersion (alpha)
        // Mean is provided but will be overwritten per observation by CaseCountLikelihood.
        RealParameter initMean = new RealParameter(new Double[] { 1.0 });
        RealParameter dispersion = new RealParameter(new Double[] { 0.5 }); // example overdispersion
        GammaPoisson nbin = new GammaPoisson(initMean, dispersion);
        nbin.initAndValidate();

        // 4) Construct the likelihood and evaluate logP
        CaseCountLikelihood llik = new CaseCountLikelihood();
        llik.initByName(
                "NeDynamics", neList,
                "caseCounts", caseData,
                "distribution", nbin
        );
        llik.initAndValidate();
        double logP = llik.calculateLogP();

        // 5) Compare with placeholder expected value; replace with manual calculation later
        double expectedLogP = -34.89995817214991;
        double tolerance = 1e-9;

        // For now, assert equality to the placeholder to wire the test. Update when true value known.
        assertEquals(expectedLogP, logP, tolerance,
                "CaseCountLikelihood logP does not match the expected placeholder value. Replace expectedLogP once the true value is known.");
    }
    @Test
    public void testTwoDemesFiveTimePointsWithSkylineLargeAlpha() throws Exception {
        // 1) Prepare case count observations: 2 demes x 5 time points
        // Times in years (relative to most recent sample)
        Double[] times = new Double[] {0.0, 0.1, 0.11, 0.24, 0.3,
                                       0.0, 0.15, 0.16, 0.25, 0.37};
        // Trait indices: first 5 are deme 0, next 5 are deme 1
        Double[] traits = new Double[] {0.0, 0.0, 0.0, 0.0, 0.0,
                                        1.0, 1.0, 1.0, 1.0, 1.0};
        // Case counts (integers represented as doubles)
        Double[] counts = new Double[] {5.0, 7.0, 6.0, 8.0, 10.0,
                                        3.0, 4.0, 5.0, 6.0, 8.0};

        CaseCountData caseData = new CaseCountData();
        caseData.initByName(
                "caseCounts", new RealParameter(counts),
                "observationTimes", new RealParameter(times),
                "traitIndices", new RealParameter(traits)
        );

        // 2) Prepare Ne skyline with native MASCOT classes using BEAST-style initByName
        // Build a tiny tree and trait set to satisfy RateShifts' expected inputs
        Sequence s1 = new Sequence();
        Sequence s2 = new Sequence();
        Sequence s3 = new Sequence();
        Sequence s4 = new Sequence();
        s1.initByName("taxon", "a1", "value", "???");
        s2.initByName("taxon", "b1", "value", "???");
        s3.initByName("taxon", "a2", "value", "???");
        s4.initByName("taxon", "b2", "value", "???");
        Alignment alignment = new Alignment();
        alignment.initByName("sequence", s1, "sequence", s2, "sequence", s3, "sequence", s4);

        TaxonSet taxa = new TaxonSet();
        taxa.initByName("alignment", alignment);

        TraitSet traitSet = new TraitSet();
        traitSet.initByName("value", "a1=Deme1,a2=Deme1,b1=Deme2,b2=Deme2", "traitname", "type", "taxa", taxa);
        
        // 4-tip ultrametric tree with root height ~0.4
        // Left clade: (a1:0.1,a2:0.1):0.3 -> total to tips = 0.4
        // Right clade: (b1:0.2,b2:0.2):0.2 -> total to tips = 0.4
        Tree tree = new TreeParser("((a1:0.1,a2:0.1):0.3,(b1:0.2,b2:0.2):0.2);");
        tree.initByName("taxonset", taxa, "trait", traitSet);

        // Instantiate via reflection to avoid hard compile dependency if the plugin isn't present yet
        Class<?> clsRateShifts = Class.forName("mascot.dynamics.RateShifts");
        Object rateShifts1 = clsRateShifts.getDeclaredConstructor().newInstance();
        Object rateShifts2 = clsRateShifts.getDeclaredConstructor().newInstance();
        clsRateShifts.getMethod("initByName", Object[].class)
                .invoke(rateShifts1, new Object[]{new Object[]{"tree", tree, "value", "0.25 0.5 0.75 1.00"}});
        clsRateShifts.getMethod("initByName", Object[].class)
                .invoke(rateShifts2, new Object[]{new Object[]{"tree", tree, "value", "0.25 0.5 0.75 1.00"}});

        Class<?> clsSkygrowth = Class.forName("mascot.parameterdynamics.Skygrowth");
        Object skygrowth1 = clsSkygrowth.getDeclaredConstructor().newInstance();
        Object skygrowth2 = clsSkygrowth.getDeclaredConstructor().newInstance();
        // Use 3 logNe values to match 3 intervals
        RealParameter logNeDeme1 = new RealParameter(new Double[]{Math.log(5.0), Math.log(7.5), Math.log(9.0), Math.log(3.0), Math.log(1.0)});
        RealParameter logNeDeme2 = new RealParameter(new Double[]{Math.log(3.0), Math.log(4.5), Math.log(6.0), Math.log(1.0), Math.log(0.5)});
        clsSkygrowth.getMethod("initByName", Object[].class)
                .invoke(skygrowth1, new Object[]{new Object[]{"logNe", logNeDeme1, "rateShifts", rateShifts1}});
        clsSkygrowth.getMethod("initByName", Object[].class)
                .invoke(skygrowth2, new Object[]{new Object[]{"logNe", logNeDeme2, "rateShifts", rateShifts2}});

        Class<?> clsNeDynList = Class.forName("mascot.util.InitializedNeDynamicsList");
        Object neList = clsNeDynList.getDeclaredConstructor().newInstance();
        // Add both Skygrowth dynamics as children
        clsNeDynList.getMethod("initByName", Object[].class)
                .invoke(neList, new Object[]{new Object[]{"neDynamics", skygrowth1, "neDynamics", skygrowth2}});

        // 3) Prepare Gamma-Poisson distribution with overdispersion (alpha)
        // Mean is provided but will be overwritten per observation by CaseCountLikelihood.
        RealParameter initMean = new RealParameter(new Double[] { 1.0 });
        RealParameter dispersion = new RealParameter(new Double[] { 2.0 }); // example overdispersion
        GammaPoisson nbin = new GammaPoisson(initMean, dispersion);
        nbin.initAndValidate();

        // 4) Construct the likelihood and evaluate logP
        CaseCountLikelihood llik = new CaseCountLikelihood();
        llik.initByName(
                "NeDynamics", neList,
                "caseCounts", caseData,
                "distribution", nbin
        );
        llik.initAndValidate();
        double logP = llik.calculateLogP();

        // 5) Compare with placeholder expected value; replace with manual calculation later
        double expectedLogP = -36.384059617326194;
        double tolerance = 1e-9;

        // For now, assert equality to the placeholder to wire the test. Update when true value known.
        assertEquals(expectedLogP, logP, tolerance,
                "CaseCountLikelihood logP does not match the expected");
    }

    // Helper: check presence of native classes to decide whether to enable the test
    static boolean mascotClassesPresent() {
        try {
            Class.forName("mascot.util.InitializedNeDynamicsList");
            Class.forName("mascot.parameterdynamics.Skygrowth");
            Class.forName("mascot.dynamics.RateShifts");
            return true;
        } catch (ClassNotFoundException e) {
            return false;
        }
    }
}
