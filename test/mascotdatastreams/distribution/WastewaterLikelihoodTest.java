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
import mascotdatastreams.dynamics.Spline;

public class WastewaterLikelihoodTest {

    @Test
    public void testSplinePrevalenceTwoDemes() throws Exception {
        Double[] times0 = new Double[] {0.0, 0.1, 0.52, 0.98, 1.47, 2.0};
        Double[] concentrations0 = new Double[] {0.0025, 0.2, 1.5, 35.7, 14.8, 0.8};
        Double[] times1 = new Double[] {0.0, 0.15, 0.56, 0.98, 1.78, 1.98};
        Double[] concentrations1 = new Double[] {0.01, 0.54, 2.3, 4.5, 1.2, 0.9};
        
        RateShifts rateShifts = buildRateShifts("0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0");
        RateShifts gridRateShifts = buildRateShifts("0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0");

        Double[] logI_deme0 = new Double[] { 0.0, 1.0, 2.0, 5.0, 10.0, 2.5, -2.0, 0.0, 1.0, -1.0, 0.0};
        Double[] logI_deme1 = new Double[] { 0.0, 1.0, 3.0, 6.0, 8.0, 3.0, -2.0, 0.0, -0.5, -1.0, 0.0};

        RealParameter uninf = new RealParameter(new Double[] {75.0});

        Spline spline0 = new Spline();
        spline0.initByName(
                "logInfected", new RealParameter(logI_deme0),
                "rateShifts", rateShifts,
                "gridRateShifts", gridRateShifts,
                "uninfectiousRate", uninf
        );
        spline0.initAndValidate();

        Spline spline1 = new Spline();
        spline1.initByName(
                "logInfected", new RealParameter(logI_deme1),
                "rateShifts", rateShifts,
                "gridRateShifts", gridRateShifts,
                "uninfectiousRate", uninf
        );
        spline1.initAndValidate();

        // Distribution with standard deviation on log scale; mean will be overwritten per obs by the likelihood
        RealParameter initMean = new RealParameter(new Double[]{1.0});
        RealParameter sd = new RealParameter(new Double[]{0.5});
        LogNormal ln = new LogNormal(initMean, sd);
        ln.initAndValidate();

        // Likelihood per deme (single-deme mode) using prevalenceSpline
        WastewaterLikelihood llikDeme0 = new WastewaterLikelihood();
        llikDeme0.initByName(
                "prevalenceSpline", spline0,
                "concentrations", new RealParameter(concentrations0),
                "concentrationTimes", new RealParameter(times0),
                "populationSize", new RealParameter(new Double[]{10000.0}),
                "distribution", ln
        );
        llikDeme0.initAndValidate();
        double logP0 = llikDeme0.calculateLogP();

        WastewaterLikelihood llikDeme1 = new WastewaterLikelihood();
        llikDeme1.initByName(
                "prevalenceSpline", spline1,
                "concentrations", new RealParameter(concentrations1),
                "concentrationTimes", new RealParameter(times1),
                "populationSize", new RealParameter(new Double[]{10000.0}),
                "distribution", ln
        );
        llikDeme1.initAndValidate();
        double logP1 = llikDeme1.calculateLogP();
        double logP = logP0 + logP1;

        assertEquals(-1577.3565860026617, logP, 1e-6, "Wastewater concentration likelihood should match expected value.");
    }
    
    @Test
    public void testSplinePrevalenceTwoDemesScaling() throws Exception {
        // Observations per deme: 
        Double[] times0 = new Double[] {0.0, 0.1, 0.52, 0.98, 1.47, 2.0};
        Double[] concentrations0 = new Double[] {0.0025, 0.2, 1.5, 35.7, 14.8, 0.8};
        Double[] times1 = new Double[] {0.0, 0.15, 0.56, 0.98, 1.78, 1.98};
        Double[] concentrations1 = new Double[] {0.01, 0.54, 2.3, 4.5, 1.2, 0.9};
        
        // Spline prevalence with fractional shifts at 0.5 and 1.0 (root height = 2.0 in helper)
        RateShifts rateShifts = buildRateShifts("0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0");
        RateShifts gridRateShifts = buildRateShifts("0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0");

        // 3 control points per deme (dimension = shifts + 1), in log-space
        Double[] logI_deme0 = new Double[] { 0.0, 1.0, 2.0, 5.0, 10.0, 2.5, -2.0, 0.0, 1.0, -1.0, 0.0};
        Double[] logI_deme1 = new Double[] { 0.0, 1.0, 3.0, 6.0, 8.0, 3.0, -2.0, 0.0, -0.5, -1.0, 0.0};

        // uninfectious rate is required by the spline but does not affect the mean directly here
        RealParameter uninf = new RealParameter(new Double[] {75.0});

        Spline spline0 = new Spline();
        spline0.initByName(
                "logInfected", new RealParameter(logI_deme0),
                "rateShifts", rateShifts,
                "gridRateShifts", gridRateShifts,
                "uninfectiousRate", uninf
        );
        spline0.initAndValidate();

        Spline spline1 = new Spline();
        spline1.initByName(
                "logInfected", new RealParameter(logI_deme1),
                "rateShifts", rateShifts,
                "gridRateShifts", gridRateShifts,
                "uninfectiousRate", uninf
        );
        spline1.initAndValidate();

        // Distribution with standard deviation on log scale; mean will be overwritten per obs by the likelihood
        RealParameter initMean = new RealParameter(new Double[]{1.0});
        RealParameter sd = new RealParameter(new Double[]{0.1});
        LogNormal ln = new LogNormal(initMean, sd);
        ln.initAndValidate();

        // Likelihood per deme (single-deme mode) using prevalenceSpline
        WastewaterLikelihood llikDeme0 = new WastewaterLikelihood();
        llikDeme0.initByName(
                "prevalenceSpline", spline0,
                "concentrations", new RealParameter(concentrations0),
                "concentrationTimes", new RealParameter(times0),
                "populationSize", new RealParameter(new Double[]{10000.0}),
                "distribution", ln,
                "scaling", new RealParameter(new Double[]{0.1})
        );
        llikDeme0.initAndValidate();
        double logP0 = llikDeme0.calculateLogP();

        WastewaterLikelihood llikDeme1 = new WastewaterLikelihood();
        llikDeme1.initByName(
                "prevalenceSpline", spline1,
                "concentrations", new RealParameter(concentrations1),
                "concentrationTimes", new RealParameter(times1),
                "populationSize", new RealParameter(new Double[]{10000.0}),
                "distribution", ln,
                "scaling", new RealParameter(new Double[]{0.05})
        );
        llikDeme1.initAndValidate();
        double logP1 = llikDeme1.calculateLogP();
        double logP = logP0 + logP1;

        assertEquals(-66701.11332701507, logP, 1e-6, "Wastewater concentration likelihood with scaling should match expected value.");
    }
    
    @Test
    public void testSplinePrevalenceTwoDemesConcentrationsOutsideTree() throws Exception {
        Double[] times0 = new Double[] {-0.1, 0.0, 0.1, 0.52, 0.98, 1.47, 2.0};
        Double[] concentrations0 = new Double[] {0.002, 0.0025, 0.2, 1.5, 35.7, 14.8, 0.8};
        Double[] times1 = new Double[] {-0.1, 0.0, 0.15, 0.56, 0.98, 1.78, 1.98};
        Double[] concentrations1 = new Double[] {0.12, 0.01, 0.54, 2.3, 4.5, 1.2, 0.9};
        
        // Spline prevalence with fractional shifts at 0.5 and 1.0 (root height = 2.0 in helper)
        RateShifts rateShifts = buildRateShifts("0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0");
        RateShifts gridRateShifts = buildRateShifts("0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0");

        // 3 control points per deme (dimension = shifts + 1), in log-space
        Double[] logI_deme0 = new Double[] { 0.0, 1.0, 2.0, 5.0, 10.0, 2.5, -2.0, 0.0, 1.0, -1.0, 0.0};
        Double[] logI_deme1 = new Double[] { 0.0, 1.0, 3.0, 6.0, 8.0, 3.0, -2.0, 0.0, -0.5, -1.0, 0.0};

        // uninfectious rate is required by the spline but does not affect the mean directly here
        RealParameter uninf = new RealParameter(new Double[] {75.0});

        Spline spline0 = new Spline();
        spline0.initByName(
                "logInfected", new RealParameter(logI_deme0),
                "rateShifts", rateShifts,
                "gridRateShifts", gridRateShifts,
                "uninfectiousRate", uninf
        );
        spline0.initAndValidate();

        Spline spline1 = new Spline();
        spline1.initByName(
                "logInfected", new RealParameter(logI_deme1),
                "rateShifts", rateShifts,
                "gridRateShifts", gridRateShifts,
                "uninfectiousRate", uninf
        );
        spline1.initAndValidate();

        // Distribution with standard deviation on log scale; mean will be overwritten per obs by the likelihood
        RealParameter initMean = new RealParameter(new Double[]{1.0});
        RealParameter sd = new RealParameter(new Double[]{0.5});
        LogNormal ln = new LogNormal(initMean, sd);
        ln.initAndValidate();

        // Likelihood per deme (single-deme mode) using prevalenceSpline
        WastewaterLikelihood llikDeme0 = new WastewaterLikelihood();
        llikDeme0.initByName(
                "prevalenceSpline", spline0,
                "concentrations", new RealParameter(concentrations0),
                "concentrationTimes", new RealParameter(times0),
                "populationSize", new RealParameter(new Double[]{10000.0}),
                "distribution", ln
        );
        llikDeme0.initAndValidate();
        double logP0 = llikDeme0.calculateLogP();

        WastewaterLikelihood llikDeme1 = new WastewaterLikelihood();
        llikDeme1.initByName(
                "prevalenceSpline", spline1,
                "concentrations", new RealParameter(concentrations1),
                "concentrationTimes", new RealParameter(times1),
                "populationSize", new RealParameter(new Double[]{10000.0}),
                "distribution", ln
        );
        llikDeme1.initAndValidate();
        double logP1 = llikDeme1.calculateLogP();
        double logP = logP0 + logP1;

        assertEquals(-1693.0659044120375, logP, 1e-6, "Wastewater concentration likelihood with observations outside tree should match expected value.");
    }



    private static RateShifts buildRateShifts(String shiftValues) throws Exception {
        RateShifts rs = new RateShifts();
        rs.initByName("value", shiftValues);
        return rs;
    }
}
