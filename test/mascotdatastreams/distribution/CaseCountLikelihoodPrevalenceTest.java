package mascotdatastreams.distribution;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

import beast.base.inference.parameter.RealParameter;
import mascotdatastreams.dynamics.PrevalenceConstant;
import mascotdatastreams.dynamics.PrevalenceDynamicsList;

public class CaseCountLikelihoodPrevalenceTest {

    @Test
    public void testConstantPrevalenceTwoDemes() throws Exception {
        // Observations: 2 demes x 5 time points
        Double[] times = new Double[] {0.0, 0.1, 0.11, 0.24, 0.3,
                                       0.0, 0.15, 0.16, 0.25, 0.37};
        Double[] traits = new Double[] {0.0, 0.0, 0.0, 0.0, 0.0,
                                        1.0, 1.0, 1.0, 1.0, 1.0};
        Double[] counts = new Double[] {5.0, 7.0, 6.0, 8.0, 10.0,
                                        3.0, 4.0, 5.0, 6.0, 8.0};

        CaseCountData caseData = new CaseCountData();
        caseData.initByName(
                "caseCounts", new RealParameter(counts),
                "observationTimes", new RealParameter(times),
                "traitIndices", new RealParameter(traits)
        );

        // Prevalence: constant log I = log(5) for both demes
        PrevalenceConstant prev1 = new PrevalenceConstant();
        PrevalenceConstant prev2 = new PrevalenceConstant();
        prev1.initByName("logPrevalence", new RealParameter(new Double[]{Math.log(5.0)}));
        prev2.initByName("logPrevalence", new RealParameter(new Double[]{Math.log(5.0)}));

        PrevalenceDynamicsList prevList = new PrevalenceDynamicsList();
        prevList.initByName("prevalence", prev1, "prevalence", prev2);

        // Distribution with dispersion alpha; mean will be overwritten per obs by the likelihood
        RealParameter initMean = new RealParameter(new Double[]{1.0});
        RealParameter alpha = new RealParameter(new Double[]{0.5});
        GammaPoisson gp = new GammaPoisson(initMean, alpha);
        gp.initAndValidate();

        // Likelihood
        CaseCountLikelihood llik = new CaseCountLikelihood();
        llik.initByName(
                "prevalence", prevList,
                "caseCounts", caseData,
                "distribution", gp
        );
        llik.initAndValidate();
        double logP = llik.calculateLogP();

        // Expected: sum of log PMF with mean = 5 for all observations
        double expected = 0.0;
        for (Double x : counts) {
            // Set mean to 5 for the distribution and compute log PMF
            gp.meanInput.setValue(new RealParameter(new Double[]{5.0}), gp);
            expected += gp.logPmf(x.intValue());
        }

        assertEquals(expected, logP, 1e-9, "Prevalence-based likelihood should match manual sum for constant prevalence.");
    }
}
