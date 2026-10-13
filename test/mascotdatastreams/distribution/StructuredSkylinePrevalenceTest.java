package mascotdatastreams.distribution;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import mascot.dynamics.RateShifts;
import mascot.parameterdynamics.NeDynamics;
import mascot.util.InitializedNeDynamicsList;
import mascotdatastreams.dynamics.Spline;
import mascotdatastreams.dynamics.SplinePrevalenceToNe;

/**
 * Tests the optional migrationIndices input, which ties several directed migration
 * routes to the same estimated rate (and indicator).
 */
public class StructuredSkylinePrevalenceTest {

    private static final String KNOTS = "0.0 0.5 1.0";
    private static final String GRID = "0.0 0.25 0.5 0.75 1.0";
    private static final double[] LOG_PREVALENCE = { 2.0, 3.0, 4.0 };

    private static RateShifts rateShifts(String values) {
        RateShifts rs = new RateShifts();
        rs.initByName("value", values);
        return rs;
    }

    private static InitializedNeDynamicsList neDynamics() {
        RateShifts knots = rateShifts(KNOTS);
        RateShifts grid = rateShifts(GRID);
        List<NeDynamics> list = new ArrayList<>();
        for (double logI : LOG_PREVALENCE) {
            Spline spline = new Spline();
            spline.initByName(
                    "logInfected", new RealParameter(new Double[] { logI, logI, logI }),
                    "rateShifts", knots,
                    "gridRateShifts", grid,
                    "uninfectiousRate", new RealParameter(new Double[] { 10.0 }),
                    "clipTransRate", false);
            SplinePrevalenceToNe ne = new SplinePrevalenceToNe();
            ne.initByName("spline", spline);
            list.add(ne);
        }
        InitializedNeDynamicsList neList = new InitializedNeDynamicsList();
        neList.initByName("neDynamics", list);
        return neList;
    }

    private static StructuredSkylinePrevalence dynamics(
            String rates, String indicators, String indices) {
        StructuredSkylinePrevalence dyn = new StructuredSkylinePrevalence();
        List<Object> args = new ArrayList<>(List.of(
                "dimension", 3,
                "NeDynamics", neDynamics(),
                "forwardsMigration", new RealParameter(rates),
                "rateShifts", rateShifts(KNOTS),
                "indicators", new BooleanParameter(indicators)));
        if (indices != null) {
            args.add("migrationIndices");
            args.add(new IntegerParameter(indices));
        }
        dyn.initByName(args.toArray());
        return dyn;
    }

    @Test
    public void tiedRatesMatchExpandedRates() {
        // Tie 1->3 to 1->2: slots (1->2, 1->3, 2->1, 2->3, 3->1, 3->2) -> 0 0 1 2 3 4.
        StructuredSkylinePrevalence tied = dynamics(
                "0.5 2.0 3.0 4.0 5.0", "true true true true true", "0 0 1 2 3 4");
        StructuredSkylinePrevalence full = dynamics(
                "0.5 0.5 2.0 3.0 4.0 5.0", "true true true true true true", null);

        assertEquals(5, tied.f_mInput.get().getDimension(), "tied rates must not be resized");
        assertEquals(5, tied.indicatorInput.get().getDimension());
        for (int i = 0; i < 2; i++) {
            assertArrayEquals(full.getBackwardsMigration(i), tied.getBackwardsMigration(i), 1e-12);
        }
    }

    @Test
    public void sharedIndicatorSwitchesAllTiedRoutes() {
        // Indicator 0 is shared by 1->2 and 1->3, so both routes are off.
        StructuredSkylinePrevalence tied = dynamics(
                "0.5 2.0 3.0 4.0 5.0", "false true true true true", "0 0 1 2 3 4");

        double[] m = tied.getBackwardsMigration(0);
        // backwards matrix entry [a*dim+b] holds forward route b->a
        assertEquals(0.0, m[1 * 3 + 0], 0.0, "forward 1->2 should be off");
        assertEquals(0.0, m[2 * 3 + 0], 0.0, "forward 1->3 should be off");
        assertEquals(4 * 2, tied.getIndicators(0).length, "4 active directed routes");
    }

    @Test
    public void wrongNumberOfIndicesThrows() {
        assertThrows(RuntimeException.class, () -> dynamics(
                "0.5 2.0 3.0 4.0 5.0", "true true true true true", "0 0 1 2 3"));
    }

    @Test
    public void unusedRateThrows() {
        assertThrows(RuntimeException.class, () -> dynamics(
                "0.5 2.0 3.0 4.0 5.0 6.0", "true true true true true true", "0 0 1 2 3 4"));
    }

    @Test
    public void indexOutOfRangeThrows() {
        assertThrows(RuntimeException.class, () -> dynamics(
                "0.5 2.0 3.0 4.0 5.0", "true true true true true", "0 0 1 2 3 5"));
    }
}
