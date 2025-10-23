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

        // Create the NotAKnotSpline first
        NotAKnotSpline spline = new NotAKnotSpline();
        spline.initByName(
                "logInfected", new RealParameter(logI),
                "rateShifts", rateShifts,
                "gridRateShifts", gridRateShifts,
                "uninfectiousRate", gamma
        );
        spline.initAndValidate();

        double tolerance = 1e-12;
        // at gridRateShift points
        assertEquals(5.0, spline.getDerivative(0.0), tolerance, "Derivative at t=0.0");
        assertEquals(5.0, spline.getDerivative(0.24), tolerance, "Derivative at t=0.24");
        assertEquals(19.413793103448313, spline.getDerivative(0.76), tolerance, "Derivative at t=0.76");
        assertEquals(-13.749999999999993, spline.getDerivative(1.68), tolerance, "Derivative at t=1.68");
        assertEquals(12.5, spline.getDerivative(2.0), tolerance, "Derivative at t=2.0");
        
        // at random times in between knots and grid points, getDerivative() does not use the grid points
        assertEquals(5.0, spline.getDerivative(0.11), tolerance, "Derivative at t=0.11");
        assertEquals(13.84482758620689, spline.getDerivative(0.46), tolerance, "Derivative at t=0.46");
        assertEquals(-47.605263157894754, spline.getDerivative(0.89), tolerance, "Derivative at t=0.89");
        assertEquals(-21.585937500000032, spline.getDerivative(1.13), tolerance, "Derivative at t=1.13");
        assertEquals(11.0, spline.getDerivative(1.98), tolerance, "Derivative at t=1.98");
    }

}

