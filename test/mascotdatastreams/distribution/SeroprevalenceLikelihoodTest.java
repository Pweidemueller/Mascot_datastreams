package mascotdatastreams.distribution;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import mascot.dynamics.RateShifts;
import mascotdatastreams.dynamics.Spline;

public class SeroprevalenceLikelihoodTest {

    @Test
    public void constantPrevalence_integral_drivesLogLikelihood() throws Exception {
        // Non-uniform grid across [0, 2.0]
        RateShifts rateShifts = buildRateShifts("0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0");
        RateShifts gridRateShifts = buildRateShifts("0.0 0.1 0.52 0.7 1.25 1.5 1.9 2.0");

        // Constant prevalence c => log(c) at all knots
        double c = 10;
        Double[] logI = new Double[rateShifts.getDimension()];
        for (int i = 0; i < logI.length; i++) logI[i] = Math.log(c);

        RealParameter gamma = new RealParameter(new Double[] { 1.0 });

        Spline spline = new Spline();
        spline.initByName(
                "logInfected", new RealParameter(logI),
                "rateShifts", rateShifts,
                "gridRateShifts", gridRateShifts,
                "uninfectiousRate", gamma
        );
        spline.initAndValidate();

        // One observation at time t
        double t = 1.25; // lies on a grid point here but logic is general
        int n = 100;
        int x = 3;

        RealParameter tested = new RealParameter(new Double[] { (double) n });
        RealParameter pos = new RealParameter(new Double[] { (double) x });
        RealParameter times = new RealParameter(new Double[] { t });

        // scaling so p = scaling * c * t ∈ (0,1)
        double scalingVal = 0.5; // p = 0.5 * 10 * 1.25 = 6.25
        RealParameter scaling = new RealParameter(new Double[] { scalingVal });

        Binomial binom = new Binomial();

        SeroprevalenceLikelihood L = new SeroprevalenceLikelihood();
        L.initByName(
                "prevalenceSpline", spline,
                "SeroPeopleTested", tested,
                "SeroPeopleSeropositive", pos,
                "SeroTimes", times,
                "distribution", binom,
                "scaling", scaling
        );

        double expectedP = Math.min(1.0 - 1e-16, Math.max(1e-16, scalingVal * c * t));
        double expectedLogP = binom.logPMFForParams(x, n, expectedP);

        double actual = L.calculateLogP();
        assertEquals(expectedLogP, actual, 1e-9);
    }

    private static RateShifts buildRateShifts(String shiftValues) throws Exception {
        RateShifts rs = new RateShifts();
        rs.initByName("value", shiftValues);
        return rs;
    }
}


