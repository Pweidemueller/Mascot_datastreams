package mascotdatastreams.dynamics;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

import beast.base.inference.parameter.RealParameter;
import mascot.dynamics.RateShifts;

public class NotAKnotSplineTest {

    private static RateShifts buildRateShifts(String shiftValues) throws Exception {
        RateShifts rs = new RateShifts();
        rs.initByName("value", shiftValues);
        return rs;
    }

    @Test
    public void testDerivative() throws Exception {
    	RateShifts rateShifts = buildRateShifts("0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0");
        RateShifts gridRateShifts = buildRateShifts("0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0");

        Double[] logI = new Double[] { 0.0, 1.0, 2.0, 5.0, 10.0, 2.5, -2.0, 0.0, 1.0, -1.0, 0.0};

        // gamma (uninfectious rate)
        RealParameter gamma = new RealParameter(new Double[] { 75.0 });

        // Create the Spline first
        Spline spline = new Spline();
        spline.initByName(
                "logInfected", new RealParameter(logI),
                "rateShifts", rateShifts,
                "gridRateShifts", gridRateShifts,
                "uninfectiousRate", gamma
        );
        spline.initAndValidate();

        double tolerance = 1e-12;
        // at random times in between knots and grid points, this returns derivative at the actual time not the nearest gridpoint
        assertEquals(5.194724946469638, spline.getDerivative(0.0), tolerance, "Derivative at t=0.0");
        assertEquals(5.158214019006581, spline.getDerivative(0.05), tolerance, "Derivative at t=0.05");
        assertEquals(11.256542599592883, spline.getDerivative(0.48), tolerance, "Derivative at t=0.48");
        assertEquals(3.680299505670264, spline.getDerivative(1.52), tolerance, "Derivative at t=1.52");
        assertEquals(-1.974174905495781, spline.getDerivative(1.82), tolerance, "Derivative at t=1.82");
        assertEquals(9.877045388458594, spline.getDerivative(2.0), tolerance, "Derivative at t=2.0");

        
        // outside the knot point interval, returns no derivative
        assertEquals(0.0, spline.getDerivative(-0.2), tolerance, "Derivative at t=-0.2");
        assertEquals(0.0, spline.getDerivative(2.02), tolerance, "Derivative at t=2.02");
        
        
    }

}

