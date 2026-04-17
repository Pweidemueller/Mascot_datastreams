package mascotdatastreams.distribution;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

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

        assertEquals(-1066.432510384534, logP, 1e-6, "Wastewater concentration likelihood should match expected value.");
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

        assertEquals(-29524.927336364974, logP, 1e-6, "Wastewater concentration likelihood with scaling should match expected value.");
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

        assertEquals(-1105.35028185745, logP, 1e-6, "Wastewater concentration likelihood with observations outside tree should match expected value.");
    }


// TO DO: outdated tests they use an old parameterisation
//     /**
//      * Verifies that the WastewaterLikelihood computes the correct logP when given the true
//      * prevalence trajectory and wastewater observations from a known simulation.
//      *
//      * Simulation: 20_2_simulation (deme index==0, Sample==0)
//      * True parameters: sigma=0.7072, scaling=135.487, populationSize=50000
//      *
//      * The spline logI knots were derived by extracting I(t) for population==I, index==0, Sample==0
//      * from 20_2_simulation.traj (step-function, backward time τ = t_mostRecentSample - t_forward).
//      * The expected logP was independently verified in Python.
//      */
//     @Test
//     public void testSimulatedDataLogPAtTrueParams() throws Exception {
//         double trueScaling = 135.48730073856368;
//         double trueSigma   = 0.7072485377202837;
//         double populationSizeVal = 50000.0;

//         // 59 wastewater observations for deme 0 (t_wastewater_frommostrecentsample, concentration)
//         Double[] concentrations0 = new Double[]{
//             0.012277638061131522, 0.0067659217904971185, 0.01682766680859667,
//             0.023353461062968334, 0.0024902972245187914, 0.004780607479303126,
//             0.017763399939565636, 0.018037214242783768, 0.02854739437202502,
//             0.02157315446103002, 0.10490035487075233, 0.1268885783865827,
//             0.1253436739860014, 0.3732107479002576, 0.35153462947026587,
//             0.18465712908523405, 0.585274678349232, 0.2888365378711012,
//             1.5344890518886924, 1.1869005445122238, 1.5102732530797724,
//             1.4427372274241812, 7.533182368676005, 3.67116000179163,
//             4.0239532204013075, 5.723909452527381, 14.241374853085853,
//             15.85307935450556, 20.44673414879529, 24.67358780566904,
//             93.6619920828741, 16.977124692157787, 16.6179804957764,
//             13.62584988938819, 36.826090733718594, 50.41464217517081,
//             19.48060580831557, 10.685576787377206, 9.682263644244058,
//             24.694122100362236, 23.139837387865946, 17.56101553508704,
//             6.428528868597167, 10.420926081204222, 8.18055007234173,
//             7.351599211186624, 9.73015425254732, 5.135356008338134,
//             5.861924027995837, 3.193314967039282, 3.1498924371939068,
//             3.4055271889738226, 0.635218428329897, 1.1827340555961614,
//             0.8725236975235521, 0.6455693276456012, 0.6734190163419967,
//             2.0754945295925284, 0.3096282982619285
//         };
//         Double[] times0 = new Double[]{
//             0.1613768541536125, 0.15863712812621525, 0.15589740209881797,
//             0.15315767607142072, 0.15041795004402345, 0.1476782240166262,
//             0.14493849798922895, 0.14219877196183167, 0.13945904593443442,
//             0.13671931990703715, 0.1339795938796399, 0.13123986785224262,
//             0.12850014182484537, 0.12576041579744812, 0.12302068977005085,
//             0.1202809637426536, 0.11754123771525633, 0.11480151168785907,
//             0.11206178566046182, 0.10932205963306454, 0.1065823336056673,
//             0.10384260757827003, 0.10110288155087277, 0.0983631555234755,
//             0.09562342949607826, 0.092883703468681, 0.09014397744128373,
//             0.08740425141388647, 0.0846645253864892, 0.08192479935909194,
//             0.0791850733316947, 0.07644534730429743, 0.07370562127690017,
//             0.0709658952495029, 0.06822616922210564, 0.0654864431947084,
//             0.06274671716731113, 0.06000699113991387, 0.057267265112516605,
//             0.05452753908511934, 0.05178781305772209, 0.04904808703032483,
//             0.04630836100292757, 0.043568634975530304, 0.04082890894813304,
//             0.03808918292073579, 0.035349456893338516, 0.03260973086594127,
//             0.029870004838544018, 0.02713027881114674, 0.024390552783749492,
//             0.021650826756352215, 0.018911100728954966, 0.016171374701557717,
//             0.01343164867416044, 0.010691922646763191, 0.007952196619365914,
//             0.005212470591968665, 0.0024727445645713886
//         };

//         // Spline logI knots derived from 20_2_simulation.traj (Sample==0, population==I, index==0).
//         // 20 evenly-spaced knots in backward time [0, 0.165]; τ=0 is the most-recent-sample time.
//         // I(τ) read as step function: last event value at forward time t_mrt - τ (t_mrt ≈ 0.1614).
//         String knotTimesStr =
//             "0.00000000 0.00868421 0.01736842 0.02605263 0.03473684 0.04342105 " +
//             "0.05210526 0.06078947 0.06947368 0.07815789 0.08684211 0.09552632 " +
//             "0.10421053 0.11289474 0.12157895 0.13026316 0.13894737 0.14763158 " +
//             "0.15631579 0.16500000";
//         Double[] logI_deme0 = new Double[]{
//             5.398163, 6.011267, 6.630683, 7.198184, 7.783641, 8.331345,
//             8.793612, 9.147933, 9.336797, 9.224342, 8.711937, 7.864036,
//             6.967909, 5.817111, 4.919981, 3.688879, 2.197225, 0.693147,
//             0.000000, 0.000000
//         };

//         RealParameter uninf = new RealParameter(new Double[]{88.5968});
//         RateShifts rateShifts = buildRateShifts(knotTimesStr);

//         Spline spline0 = new Spline();
//         spline0.initByName(
//                 "logInfected", new RealParameter(logI_deme0),
//                 "rateShifts", rateShifts,
//                 "gridRateShifts", rateShifts,
//                 "uninfectiousRate", uninf
//         );
//         spline0.initAndValidate();

//         RealParameter sdParam = new RealParameter(new Double[]{trueSigma});
//         RealParameter initMean = new RealParameter(new Double[]{1.0});
//         LogNormal ln = new LogNormal(initMean, sdParam);
//         ln.initAndValidate();

//         WastewaterLikelihood wl = new WastewaterLikelihood();
//         wl.initByName(
//                 "prevalenceSpline", spline0,
//                 "concentrations", new RealParameter(concentrations0),
//                 "concentrationTimes", new RealParameter(times0),
//                 "populationSize", new RealParameter(new Double[]{populationSizeVal}),
//                 "distribution", ln,
//                 "scaling", new RealParameter(new Double[]{trueScaling})
//         );
//         wl.initAndValidate();

//         double logP = wl.calculateLogP();
//         assertEquals(-69.406190, logP, 1e-4,
//                 "logP at true parameters should match expected value");
//     }

//     /**
//      * Profiles the likelihood over sigma while holding all other parameters fixed at their
//      * true simulation values (true prevalence from traj, true scaling=135.487, true N=50000).
//      *
//      * The MLE sigma is expected to be substantially below the true sigma (0.707), demonstrating
//      * the downward bias seen in the simulation study. This test documents that the bias persists
//      * even with the exact true prevalence trajectory — pointing to the likelihood formulation
//      * (e.g. the epsilon=1e-2 floor on scaledMean, or the LogNormal parameterisation) as the source.
//      */
//     @Test
//     public void testSimulatedDataSigmaProfileShowsBias() throws Exception {
//         double trueScaling = 135.48730073856368;
//         double trueSigma   = 0.7072485377202837;
//         double populationSizeVal = 50000.0;

//         Double[] concentrations0 = new Double[]{
//             0.012277638061131522, 0.0067659217904971185, 0.01682766680859667,
//             0.023353461062968334, 0.0024902972245187914, 0.004780607479303126,
//             0.017763399939565636, 0.018037214242783768, 0.02854739437202502,
//             0.02157315446103002, 0.10490035487075233, 0.1268885783865827,
//             0.1253436739860014, 0.3732107479002576, 0.35153462947026587,
//             0.18465712908523405, 0.585274678349232, 0.2888365378711012,
//             1.5344890518886924, 1.1869005445122238, 1.5102732530797724,
//             1.4427372274241812, 7.533182368676005, 3.67116000179163,
//             4.0239532204013075, 5.723909452527381, 14.241374853085853,
//             15.85307935450556, 20.44673414879529, 24.67358780566904,
//             93.6619920828741, 16.977124692157787, 16.6179804957764,
//             13.62584988938819, 36.826090733718594, 50.41464217517081,
//             19.48060580831557, 10.685576787377206, 9.682263644244058,
//             24.694122100362236, 23.139837387865946, 17.56101553508704,
//             6.428528868597167, 10.420926081204222, 8.18055007234173,
//             7.351599211186624, 9.73015425254732, 5.135356008338134,
//             5.861924027995837, 3.193314967039282, 3.1498924371939068,
//             3.4055271889738226, 0.635218428329897, 1.1827340555961614,
//             0.8725236975235521, 0.6455693276456012, 0.6734190163419967,
//             2.0754945295925284, 0.3096282982619285
//         };
//         Double[] times0 = new Double[]{
//             0.1613768541536125, 0.15863712812621525, 0.15589740209881797,
//             0.15315767607142072, 0.15041795004402345, 0.1476782240166262,
//             0.14493849798922895, 0.14219877196183167, 0.13945904593443442,
//             0.13671931990703715, 0.1339795938796399, 0.13123986785224262,
//             0.12850014182484537, 0.12576041579744812, 0.12302068977005085,
//             0.1202809637426536, 0.11754123771525633, 0.11480151168785907,
//             0.11206178566046182, 0.10932205963306454, 0.1065823336056673,
//             0.10384260757827003, 0.10110288155087277, 0.0983631555234755,
//             0.09562342949607826, 0.092883703468681, 0.09014397744128373,
//             0.08740425141388647, 0.0846645253864892, 0.08192479935909194,
//             0.0791850733316947, 0.07644534730429743, 0.07370562127690017,
//             0.0709658952495029, 0.06822616922210564, 0.0654864431947084,
//             0.06274671716731113, 0.06000699113991387, 0.057267265112516605,
//             0.05452753908511934, 0.05178781305772209, 0.04904808703032483,
//             0.04630836100292757, 0.043568634975530304, 0.04082890894813304,
//             0.03808918292073579, 0.035349456893338516, 0.03260973086594127,
//             0.029870004838544018, 0.02713027881114674, 0.024390552783749492,
//             0.021650826756352215, 0.018911100728954966, 0.016171374701557717,
//             0.01343164867416044, 0.010691922646763191, 0.007952196619365914,
//             0.005212470591968665, 0.0024727445645713886
//         };
//         String knotTimesStr =
//             "0.00000000 0.00868421 0.01736842 0.02605263 0.03473684 0.04342105 " +
//             "0.05210526 0.06078947 0.06947368 0.07815789 0.08684211 0.09552632 " +
//             "0.10421053 0.11289474 0.12157895 0.13026316 0.13894737 0.14763158 " +
//             "0.15631579 0.16500000";
//         Double[] logI_deme0 = new Double[]{
//             5.398163, 6.011267, 6.630683, 7.198184, 7.783641, 8.331345,
//             8.793612, 9.147933, 9.336797, 9.224342, 8.711937, 7.864036,
//             6.967909, 5.817111, 4.919981, 3.688879, 2.197225, 0.693147,
//             0.000000, 0.000000
//         };

//         RealParameter uninf = new RealParameter(new Double[]{88.5968});
//         RateShifts rateShifts = buildRateShifts(knotTimesStr);

//         Spline spline0 = new Spline();
//         spline0.initByName(
//                 "logInfected", new RealParameter(logI_deme0),
//                 "rateShifts", rateShifts,
//                 "gridRateShifts", rateShifts,
//                 "uninfectiousRate", uninf
//         );
//         spline0.initAndValidate();

//         // sdParam is shared with LogNormal; we update it in the scan loop
//         RealParameter sdParam = new RealParameter(new Double[]{trueSigma});
//         RealParameter initMean = new RealParameter(new Double[]{1.0});
//         LogNormal ln = new LogNormal(initMean, sdParam);
//         ln.initAndValidate();

//         WastewaterLikelihood wl = new WastewaterLikelihood();
//         wl.initByName(
//                 "prevalenceSpline", spline0,
//                 "concentrations", new RealParameter(concentrations0),
//                 "concentrationTimes", new RealParameter(times0),
//                 "populationSize", new RealParameter(new Double[]{populationSizeVal}),
//                 "distribution", ln,
//                 "scaling", new RealParameter(new Double[]{trueScaling})
//         );
//         wl.initAndValidate();

//         // Profile likelihood over sigma in [0.05, 2.0]
//         double maxLogP = Double.NEGATIVE_INFINITY;
//         double mleSigma = -1.0;
//         int nSteps = 200;
//         for (int i = 0; i < nSteps; i++) {
//             double sigma = 0.05 + i * 0.01;
//             sdParam.setValue(0, sigma);
//             double logP = wl.calculateLogP();
//             if (logP > maxLogP) {
//                 maxLogP = logP;
//                 mleSigma = sigma;
//             }
//         }

//         System.out.printf("True sigma: %.4f%n", trueSigma);
//         System.out.printf("MLE sigma (profile at true scaling): %.3f  logP=%.4f%n", mleSigma, maxLogP);

//         // MLE sigma should be well below the true value (demonstrating the downward bias)
//         assertTrue(mleSigma < trueSigma - 0.05,
//                 String.format("MLE sigma (%.3f) should be substantially below true sigma (%.4f)", mleSigma, trueSigma));
//         // Expected MLE is ~0.58 based on Python verification
//         assertEquals(0.58, mleSigma, 0.15,
//                 "MLE sigma expected around 0.58 (biased low relative to true 0.707)");
//     }

    private static RateShifts buildRateShifts(String shiftValues) throws Exception {
        RateShifts rs = new RateShifts();
        rs.initByName("value", shiftValues);
        return rs;
    }
}
