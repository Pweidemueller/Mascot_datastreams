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
        RateShifts rateShifts = buildRateShifts("0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0");
        RateShifts gridRateShifts = buildRateShifts("0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0");

        Double[] logI = new Double[] { 0.0, 1.0, 2.0, 5.0, 10.0, 2.5, -2.0, 0.0, 1.0, 2.0, 0.0};

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
        double t = 1.25;
        int n = 100;
        int x = 3;

        RealParameter tested = new RealParameter(new Double[] { (double) n });
        RealParameter pos = new RealParameter(new Double[] { (double) x });
        RealParameter times = new RealParameter(new Double[] { t });
        RealParameter populationSize = new RealParameter(new Double[] { 10000.0 });

        // scaling so p = scaling * cumulativeIncidence / populationSize ∈ (0,1)
        // For constant prevalence c and transmission rate = 1, cumulative incidence from 0 to t is c * t
        // So p = scaling * (c * t) / populationSize
        double scalingVal = 0.5;
        RealParameter scaling = new RealParameter(new Double[] { scalingVal });

        Binomial binom = new Binomial();

        SeroprevalenceLikelihood L = new SeroprevalenceLikelihood();
        L.initByName(
                "prevalenceSpline", spline,
                "seroPeopleTested", tested,
                "seroPeopleSeropositive", pos,
                "seroTimes", times,
                "populationSize", populationSize,
                "distribution", binom,
                "scaling", scaling
        );

        double actual = L.calculateLogP();
        assertEquals(-16.971882317175133, actual, 1e-9);
    }

    private static RateShifts buildRateShifts(String shiftValues) throws Exception {
        RateShifts rs = new RateShifts();
        rs.initByName("value", shiftValues);
        return rs;
    }
}


