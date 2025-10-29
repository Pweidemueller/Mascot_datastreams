package mascotdatastreams.distribution;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

import beast.base.inference.parameter.RealParameter;
import mascot.dynamics.RateShifts;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import mascotdatastreams.dynamics.NotAKnotSpline;

public class CaseCountLikelihoodTest {

    @Test
    public void testSplinePrevalenceTwoDemes() throws Exception {
         Double[] times0 = new Double[] {0.0, 0.1, 0.52, 0.98, 1.47, 2.0};
         Double[] counts0 = new Double[] {5.0, 7.0, 6.0, 8.0, 10.0, 0.0};
         Double[] times1 = new Double[] {0.0, 0.15, 0.56, 0.98, 1.78, 1.98};
         Double[] counts1 = new Double[] {3.0, 4.0, 5.0, 6.0, 8.0, 2.0};
        
        RateShifts rateShifts = buildRateShifts("0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0");
        RateShifts gridRateShifts = buildRateShifts("0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0");

        Double[] logI_deme0 = new Double[] { 0.0, 1.0, 2.0, 5.0, 10.0, 2.5, -2.0, 0.0, 1.0, -1.0, 0.0};
        Double[] logI_deme1 = new Double[] { 0.0, 1.0, 3.0, 6.0, 8.0, 3.0, -2.0, 0.0, -0.5, -1.0, 0.0};

        RealParameter uninf = new RealParameter(new Double[] {75.0});

        NotAKnotSpline spline0 = new NotAKnotSpline();
        spline0.initByName(
                "logInfected", new RealParameter(logI_deme0),
                "rateShifts", rateShifts,
                "gridRateShifts", gridRateShifts,
                "uninfectiousRate", uninf
        );
        spline0.initAndValidate();

        NotAKnotSpline spline1 = new NotAKnotSpline();
        spline1.initByName(
                "logInfected", new RealParameter(logI_deme1),
                "rateShifts", rateShifts,
                "gridRateShifts", gridRateShifts,
                "uninfectiousRate", uninf
        );
        spline1.initAndValidate();

        // Distribution with dispersion alpha; mean will be overwritten per obs by the likelihood
        RealParameter initMean = new RealParameter(new Double[]{1.0});
        RealParameter alpha = new RealParameter(new Double[]{0.5});
        GammaPoisson gp = new GammaPoisson(initMean, alpha);
        gp.initAndValidate();

        // Likelihood per deme (single-deme mode) using prevalenceSpline
        CaseCountLikelihood llikDeme0 = new CaseCountLikelihood();
        llikDeme0.initByName(
                "prevalenceSpline", spline0,
                "caseCounts", new RealParameter(counts0),
                "caseTimes", new RealParameter(times0),
                "distribution", gp
        );
        llikDeme0.initAndValidate();
        double logP0 = llikDeme0.calculateLogP();

        CaseCountLikelihood llikDeme1 = new CaseCountLikelihood();
        llikDeme1.initByName(
                "prevalenceSpline", spline1,
                "caseCounts", new RealParameter(counts1),
                "caseTimes", new RealParameter(times1),
                "distribution", gp
        );
        llikDeme1.initAndValidate();
        double logP1 = llikDeme1.calculateLogP();
        double logP = logP0 + logP1;

        assertEquals(-53.950267384878, logP, 1e-9, "Prevalence-based likelihood should match manual sum for Mascot Skygrowth prevalence.");
    }
    
    @Test
    public void testSplinePrevalenceTwoDemesScaling() throws Exception {
        // Observations per deme: 5 time points each
        Double[] times0 = new Double[] {0.0, 0.1, 0.52, 0.98, 1.47, 2.0};
        Double[] counts0 = new Double[] {5.0, 7.0, 6.0, 8.0, 10.0, 0.0};
        Double[] times1 = new Double[] {0.0, 0.15, 0.56, 0.98, 1.78, 1.98};
        Double[] counts1 = new Double[] {3.0, 4.0, 5.0, 6.0, 8.0, 2.0};
        
        // Spline prevalence with fractional shifts at 0.5 and 1.0 (root height = 2.0 in helper)
        RateShifts rateShifts = buildRateShifts("0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0");
        RateShifts gridRateShifts = buildRateShifts("0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0");

        // 3 control points per deme (dimension = shifts + 1), in log-space
        Double[] logI_deme0 = new Double[] { 0.0, 1.0, 2.0, 5.0, 10.0, 2.5, -2.0, 0.0, 1.0, -1.0, 0.0};
        Double[] logI_deme1 = new Double[] { 0.0, 1.0, 3.0, 6.0, 8.0, 3.0, -2.0, 0.0, -0.5, -1.0, 0.0};

        // uninfectious rate is required by the spline but does not affect the mean directly here
        RealParameter uninf = new RealParameter(new Double[] {75.0});

        NotAKnotSpline spline0 = new NotAKnotSpline();
        spline0.initByName(
                "logInfected", new RealParameter(logI_deme0),
                "rateShifts", rateShifts,
                "gridRateShifts", gridRateShifts,
                "uninfectiousRate", uninf
        );
        spline0.initAndValidate();

        NotAKnotSpline spline1 = new NotAKnotSpline();
        spline1.initByName(
                "logInfected", new RealParameter(logI_deme1),
                "rateShifts", rateShifts,
                "gridRateShifts", gridRateShifts,
                "uninfectiousRate", uninf
        );
        spline1.initAndValidate();

        // Distribution with dispersion alpha; mean will be overwritten per obs by the likelihood
        RealParameter initMean = new RealParameter(new Double[]{1.0});
        RealParameter alpha = new RealParameter(new Double[]{0.1});
        GammaPoisson gp = new GammaPoisson(initMean, alpha);
        gp.initAndValidate();

        // Likelihood per deme (single-deme mode) using prevalenceSpline
        CaseCountLikelihood llikDeme0 = new CaseCountLikelihood();
        llikDeme0.initByName(
                "prevalenceSpline", spline0,
                "caseCounts", new RealParameter(counts0),
                "caseTimes", new RealParameter(times0),
                "distribution", gp,
                "scaling", new RealParameter(new Double[]{0.1})
        );
        llikDeme0.initAndValidate();
        double logP0 = llikDeme0.calculateLogP();

        CaseCountLikelihood llikDeme1 = new CaseCountLikelihood();
        llikDeme1.initByName(
                "prevalenceSpline", spline1,
                "caseCounts", new RealParameter(counts1),
                "caseTimes", new RealParameter(times1),
                "distribution", gp,
                "scaling", new RealParameter(new Double[]{0.05})
        );
        llikDeme1.initAndValidate();
        double logP1 = llikDeme1.calculateLogP();
        double logP = logP0 + logP1;

        assertEquals(-156.616656999962, logP, 1e-9, "Prevalence-based likelihood should match manual sum for Mascot Skygrowth prevalence.");
    }
    
    @Test
    public void testSplinePrevalenceTwoDemesCaseCountsOutsideTree() throws Exception {
        // Observations per deme: 5 time points each
        Double[] times0 = new Double[] {-0.1, 0.0, 0.1, 0.52, 0.98, 1.47, 2.0};
        Double[] counts0 = new Double[] {2.0, 100.0, 5480.0, 1500.0, 560.0, 34.0, 0.0};
        Double[] times1 = new Double[] {-0.1, 0.0, 0.15, 0.56, 0.98, 1.78, 1.98};
        Double[] counts1 = new Double[] {1.0, 124.0, 178.0, 1000.0, 1487.0, 246.0, 2.0};
        
        // Spline prevalence with fractional shifts at 0.5 and 1.0 (root height = 2.0 in helper)
        RateShifts rateShifts = buildRateShifts("0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0");
        RateShifts gridRateShifts = buildRateShifts("0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0");

        // 3 control points per deme (dimension = shifts + 1), in log-space
        Double[] logI_deme0 = new Double[] { 0.0, 1.0, 2.0, 5.0, 10.0, 2.5, -2.0, 0.0, 1.0, -1.0, 0.0};
        Double[] logI_deme1 = new Double[] { 0.0, 1.0, 3.0, 6.0, 8.0, 3.0, -2.0, 0.0, -0.5, -1.0, 0.0};

        // uninfectious rate is required by the spline but does not affect the mean directly here
        RealParameter uninf = new RealParameter(new Double[] {75.0});

        NotAKnotSpline spline0 = new NotAKnotSpline();
        spline0.initByName(
                "logInfected", new RealParameter(logI_deme0),
                "rateShifts", rateShifts,
                "gridRateShifts", gridRateShifts,
                "uninfectiousRate", uninf
        );
        spline0.initAndValidate();

        NotAKnotSpline spline1 = new NotAKnotSpline();
        spline1.initByName(
                "logInfected", new RealParameter(logI_deme1),
                "rateShifts", rateShifts,
                "gridRateShifts", gridRateShifts,
                "uninfectiousRate", uninf
        );
        spline1.initAndValidate();

        // Distribution with dispersion alpha; mean will be overwritten per obs by the likelihood
        RealParameter initMean = new RealParameter(new Double[]{1.0});
        RealParameter alpha = new RealParameter(new Double[]{0.5});
        GammaPoisson gp = new GammaPoisson(initMean, alpha);
        gp.initAndValidate();

        // Likelihood per deme (single-deme mode) using prevalenceSpline
        CaseCountLikelihood llikDeme0 = new CaseCountLikelihood();
        llikDeme0.initByName(
                "prevalenceSpline", spline0,
                "caseCounts", new RealParameter(counts0),
                "caseTimes", new RealParameter(times0),
                "distribution", gp
        );
        llikDeme0.initAndValidate();
        double logP0 = llikDeme0.calculateLogP();

        CaseCountLikelihood llikDeme1 = new CaseCountLikelihood();
        llikDeme1.initByName(
                "prevalenceSpline", spline1,
                "caseCounts", new RealParameter(counts1),
                "caseTimes", new RealParameter(times1),
                "distribution", gp
        );
        llikDeme1.initAndValidate();
        double logP1 = llikDeme1.calculateLogP();
        double logP = logP0 + logP1;

        assertEquals(-5543.832990969199, logP, 1e-9, "Prevalence-based likelihood should match manual sum for Mascot Skygrowth prevalence.");
    }



    private static RateShifts buildRateShifts(String shiftValues) throws Exception {
        RateShifts rs = new RateShifts();
        rs.initByName("value", shiftValues);
        return rs;
    }
}
