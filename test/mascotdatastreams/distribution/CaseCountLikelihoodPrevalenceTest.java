package mascotdatastreams.distribution;

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

public class CaseCountLikelihoodPrevalenceTest {

    @Test
    public void testSkygrowthPrevalenceTwoDemes() throws Exception {
        // Observations per deme: 5 time points each
        Double[] times0 = new Double[] {0.0, 0.1, 0.11, 0.24, 0.3};
        Double[] counts0 = new Double[] {5.0, 7.0, 6.0, 8.0, 10.0};
        Double[] times1 = new Double[] {0.0, 0.15, 0.16, 0.25, 0.37};
        Double[] counts1 = new Double[] {3.0, 4.0, 5.0, 6.0, 8.0};

        // Prevalence: Mascot Skygrowth with fractional shifts at 0.5 and 1.0 (root height = 2.0 in helper)
        Object rateShifts = buildRateShifts("0.5 1.0");
        // 3 control points per deme (dimension = shifts + 1), in log-space
        Double[] logI_deme0 = new Double[] { Math.log(5.0), Math.log(10.0), Math.log(3.0) };
        Double[] logI_deme1 = new Double[] { Math.log(3.0), Math.log(4.0), Math.log(2.0) };

        Skygrowth skygrowth0 = new Skygrowth();
        skygrowth0.initByName(
                "logNe", new RealParameter(logI_deme0),
                "rateShifts", rateShifts
        );
        skygrowth0.initAndValidate();
        Skygrowth skygrowth1 = new Skygrowth();
        skygrowth1.initByName(
                "logNe", new RealParameter(logI_deme1),
                "rateShifts", rateShifts
        );
        skygrowth1.initAndValidate();

        // Distribution with dispersion alpha; mean will be overwritten per obs by the likelihood
        RealParameter initMean = new RealParameter(new Double[]{1.0});
        RealParameter alpha = new RealParameter(new Double[]{0.5});
        GammaPoisson gp = new GammaPoisson(initMean, alpha);
        gp.initAndValidate();

        // Likelihood per deme (single-deme mode)
        CaseCountLikelihood llikDeme0 = new CaseCountLikelihood();
        llikDeme0.initByName(
                "prevalence", skygrowth0,
                "caseCounts", new RealParameter(counts0),
                "caseTimes", new RealParameter(times0),
                "distribution", gp
        );
        llikDeme0.initAndValidate();
        double logP0 = llikDeme0.calculateLogP();

        CaseCountLikelihood llikDeme1 = new CaseCountLikelihood();
        llikDeme1.initByName(
                "prevalence", skygrowth1,
                "caseCounts", new RealParameter(counts1),
                "caseTimes", new RealParameter(times1),
                "distribution", gp
        );
        llikDeme1.initAndValidate();
        double logP1 = llikDeme1.calculateLogP();
        double logP = logP0 + logP1;

        // Expected: sum of log PMF with mean = I(t) per observation using the appropriate deme's prevalence
        double expected = 0.0;
        for (int i = 0; i < counts0.length; i++) {
            double meanI = skygrowth0.getNeTime(times0[i]);
            gp.meanInput.setValue(new RealParameter(new Double[]{meanI}), gp);
            expected += gp.logPmf(counts0[i].intValue());
        }
        for (int i = 0; i < counts1.length; i++) {
            double meanI = skygrowth1.getNeTime(times1[i]);
            gp.meanInput.setValue(new RealParameter(new Double[]{meanI}), gp);
            expected += gp.logPmf(counts1[i].intValue());
        }

        assertEquals(expected, logP, 1e-9, "Prevalence-based likelihood should match manual sum for Mascot Skygrowth prevalence.");
    }

    @Test
    public void testSkygrowthPrevalenceTwoDemesWithScaling() throws Exception {
        // Observations per deme: 5 time points each
        Double[] times0 = new Double[] {0.0, 0.1, 0.11, 0.24, 0.3};
        Double[] counts0 = new Double[] {5.0, 7.0, 6.0, 8.0, 10.0};
        Double[] times1 = new Double[] {0.0, 0.15, 0.16, 0.25, 0.37};
        Double[] counts1 = new Double[] {3.0, 4.0, 5.0, 6.0, 8.0};

        Object rateShifts = buildRateShifts("0.5 1.0");
        Double[] logI_deme0 = new Double[] { Math.log(5.0), Math.log(10.0), Math.log(3.0) };
        Double[] logI_deme1 = new Double[] { Math.log(3.0), Math.log(4.0), Math.log(2.0) };
        
        Skygrowth skygrowth0 = new Skygrowth();
        skygrowth0.initByName(
                "logNe", new RealParameter(logI_deme0),
                "rateShifts", rateShifts
        );
        skygrowth0.initAndValidate();
        Skygrowth skygrowth1 = new Skygrowth();
        skygrowth1.initByName(
                "logNe", new RealParameter(logI_deme1),
                "rateShifts", rateShifts
        );
        skygrowth1.initAndValidate();

        RealParameter initMean = new RealParameter(new Double[]{1.0});
        RealParameter alpha = new RealParameter(new Double[]{0.5});
        GammaPoisson gp = new GammaPoisson(initMean, alpha);
        gp.initAndValidate();

        // scaling parameter
        RealParameter scaling = new RealParameter(new Double[]{0.1});

        CaseCountLikelihood llikDeme0 = new CaseCountLikelihood();
        llikDeme0.initByName(
                "prevalence", skygrowth0,
                "caseCounts", new RealParameter(counts0),
                "caseTimes", new RealParameter(times0),
                "distribution", gp,
                "scaling", scaling
        );
        llikDeme0.initAndValidate();
        double logP0 = llikDeme0.calculateLogP();

        CaseCountLikelihood llikDeme1 = new CaseCountLikelihood();
        llikDeme1.initByName(
                "prevalence", skygrowth1,
                "caseCounts", new RealParameter(counts1),
                "caseTimes", new RealParameter(times1),
                "distribution", gp,
                "scaling", scaling
        );
        llikDeme1.initAndValidate();
        double logP1 = llikDeme1.calculateLogP();
        double logP = logP0 + logP1;

        // Expected: sum with mean = I(t) * 0.1 per observation using the appropriate deme's prevalence
        double expected = 0.0;
        for (int i = 0; i < counts0.length; i++) {
            double meanI = skygrowth0.getNeTime(times0[i]) * 0.1;
            gp.meanInput.setValue(new RealParameter(new Double[]{meanI}), gp);
            expected += gp.logPmf(counts0[i].intValue());
        }
        for (int i = 0; i < counts1.length; i++) {
            double meanI = skygrowth1.getNeTime(times1[i]) * 0.1;
            gp.meanInput.setValue(new RealParameter(new Double[]{meanI}), gp);
            expected += gp.logPmf(counts1[i].intValue());
        }

        assertEquals(expected, logP, 1e-9, "Prevalence-based likelihood with scaling should match manual sum for Mascot Skygrowth prevalence.");
    }

    // Build a minimal tree and RateShifts with fractional breakpoints resolved against root height
    private static Object buildRateShifts(String shiftValues) throws Exception {
        // Minimal BEAST setup
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

        // 4-tip ultrametric tree; root height arbitrary (2.0)
        Tree tree = new TreeParser("((a1:1.0,a2:1.0):1.0,(b1:1.0,b2:1.0):1.0);");
        tree.initByName("taxonset", taxa, "trait", traitSet);

        Class<?> clsRateShifts = Class.forName("mascot.dynamics.RateShifts");
        Object rs = clsRateShifts.getDeclaredConstructor().newInstance();
        clsRateShifts.getMethod("initByName", Object[].class)
                .invoke(rs, new Object[] { new Object[] { "tree", tree, "value", shiftValues } });
        return rs;
    }
}
