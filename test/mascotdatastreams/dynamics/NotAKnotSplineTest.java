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
        assertEquals(6.453884388807061, spline.getDerivative(0.0, 0), tolerance, "Prevalence at t=0.0");
        assertEquals(6.453884388807062, spline.getDerivative(0.4, 2), tolerance, "Prevalence at t=0.4");
        assertEquals(-6.0995029455080925, spline.getDerivative(0.8, 4), tolerance, "Prevalence at t=0.8");
        assertEquals(-1.8469256259204478, spline.getDerivative(1.2, 6), tolerance, "Prevalence at t=1.2000000000000002");
        assertEquals(-4.7574558173784975, spline.getDerivative(1.6, 8), tolerance, "Prevalence at t=1.6");
        assertEquals(25.242544182621508, spline.getDerivative(2.0, 9), tolerance, "Prevalence at t=2.0");
        
        // at random times in between knots and grid points
        assertEquals(4.714675188696615, spline.getDerivative(0.11, 0), tolerance, "Prevalence at t=0.11");
        assertEquals(9.481008376288646, spline.getDerivative(0.46, 2), tolerance, "Prevalence at t=0.46");
        assertEquals(-41.93192827227541, spline.getDerivative(0.89, 4), tolerance, "Prevalence at t=0.89");
        assertEquals(-16.350521792157604, spline.getDerivative(1.13, 5), tolerance, "Prevalence at t=1.13");
        assertEquals(20.110919090574374, spline.getDerivative(1.98, 9), tolerance, "Prevalence at t=1.98");
    }

//    @Test
//    public void testDerivativeMultipleKnotsLinear() throws Exception {
//        RateShifts rateShifts = buildRateShifts("0.0 0.5 1.0 1.5 2.0");
//        RateShifts gridRateShifts = buildRateShifts("0.0 0.5 1.0 1.5 2.0");
//
//        double m = -0.8;
//        double b = 0.3;
//
//        Double[] logI = new Double[] {
//                b,
//                m * 0.5 + b,
//                m * 1.0 + b,
//                m * 1.5 + b,
//                m * 2.0 + b
//        };
//        RealParameter gamma = new RealParameter(new Double[] { 5.0 });
//
//        NotAKnotSpline spline = new NotAKnotSpline();
//        spline.initByName(
//                "logInfected", new RealParameter(logI),
//                "rateShifts", rateShifts,
//                "gridRateShifts", gridRateShifts,
//                "uninfectiousRate", gamma
//        );
//        spline.initAndValidate();
//
//        double tol = 1e-9;
//        assertEquals(m, spline.getDerivativeAtGridPoint(0.0), tol);
//        assertEquals(m, spline.getDerivativeAtGridPoint(0.5), tol);
//        assertEquals(m, spline.getDerivativeAtGridPoint(1.0), tol);
//        assertEquals(m, spline.getDerivativeAtGridPoint(1.5), tol);
//        assertEquals(m, spline.getDerivativeAtGridPoint(2.0), tol);
//        assertEquals(m, spline.getDerivativeAtGridPoint(0.25), tol);
//        assertEquals(m, spline.getDerivativeAtGridPoint(0.75), tol);
//        assertEquals(m, spline.getDerivativeAtGridPoint(1.25), tol);
//        assertEquals(m, spline.getDerivativeAtGridPoint(1.75), tol);
//    }
}

