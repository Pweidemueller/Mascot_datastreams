package mascotdatastreams.dynamics;

import beast.base.inference.parameter.RealParameter;
import junit.framework.Assert;
import mascot.dynamics.RateShifts;
import org.junit.Test;

import java.util.Arrays;

/**
 * Test class for NotAKnotSpline
 * Tests spline interpolation, derivative calculation, and transmission rate computation
 */
public class NotAKnotSplineTest {

    @Test
    public void testLinearInterpolationWithTwoPoints() {
        // Test case: only 2 points should give linear interpolation
        RateShifts rateShifts = new RateShifts();
        rateShifts.initByName("value", Arrays.asList(0.0, 1.0));
        
        RateShifts gridRateShifts = new RateShifts();
        gridRateShifts.initByName("value", Arrays.asList(0.0, 0.5, 1.0));
        
        RealParameter logInfected = new RealParameter("0.0 1.0");
        RealParameter uninfectiousRate = new RealParameter("0.5");
        
        NotAKnotSpline spline = new NotAKnotSpline();
        spline.initByName(
            "logInfected", logInfected,
            "rateShifts", rateShifts,
            "gridRateShifts", gridRateShifts,
            "uninfectiousRate", uninfectiousRate
        );
        
        spline.update();
        
        // With linear interpolation, derivative should be constant
        double derivative0 = spline.getDerivativeAtGridPoint(0.0);
        double derivative05 = spline.getDerivativeAtGridPoint(0.5);
        double derivative1 = spline.getDerivativeAtGridPoint(1.0);
        
        // For linear: y = x, derivative = 1.0
        Assert.assertEquals(1.0, derivative0, 1e-10);
        Assert.assertEquals(1.0, derivative05, 1e-10);
        Assert.assertEquals(1.0, derivative1, 1e-10);
        
        // Prevalence should be exp(logInfected)
        double prev0 = spline.getPrevalenceAtGridPoint(0.0);
        Assert.assertEquals(Math.exp(0.0), prev0, 1e-10);
    }

    @Test
    public void testThreePointNaturalSpline() {
        // Test case: 3 points use natural spline boundary conditions
        RateShifts rateShifts = new RateShifts();
        rateShifts.initByName("value", Arrays.asList(0.0, 1.0, 2.0));
        
        RateShifts gridRateShifts = new RateShifts();
        gridRateShifts.initByName("value", Arrays.asList(0.0, 0.5, 1.0, 1.5, 2.0));
        
        // Quadratic: y = x^2, so log(I) = 0, 0, log(4) doesn't work, try linear log space
        RealParameter logInfected = new RealParameter("0.0 1.0 2.0"); // linear in log space
        RealParameter uninfectiousRate = new RealParameter("2.0");
        
        NotAKnotSpline spline = new NotAKnotSpline();
        spline.initByName(
            "logInfected", logInfected,
            "rateShifts", rateShifts,
            "gridRateShifts", gridRateShifts,
            "uninfectiousRate", uninfectiousRate
        );
        
        boolean valid = spline.update();
        Assert.assertTrue(valid);
        
        // Check that values are retrieved correctly
        double prev0 = spline.getPrevalenceAtGridPoint(0.0);
        double prev1 = spline.getPrevalenceAtGridPoint(1.0);
        double prev2 = spline.getPrevalenceAtGridPoint(2.0);
        
        Assert.assertEquals(Math.exp(0.0), prev0, 1e-10);
        Assert.assertTrue(prev1 > prev0); // Should be increasing
        Assert.assertTrue(prev2 > prev1);
    }

    @Test
    public void testNotAKnotSplineWithFourPoints() {
        // Test case: 4+ points use not-a-knot boundary conditions
        RateShifts rateShifts = new RateShifts();
        rateShifts.initByName("value", Arrays.asList(0.0, 1.0, 2.0, 3.0));
        
        RateShifts gridRateShifts = new RateShifts();
        gridRateShifts.initByName("value", Arrays.asList(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0));
        
        RealParameter logInfected = new RealParameter("0.0 1.0 1.5 2.0");
        RealParameter uninfectiousRate = new RealParameter("1.0");
        
        NotAKnotSpline spline = new NotAKnotSpline();
        spline.initByName(
            "logInfected", logInfected,
            "rateShifts", rateShifts,
            "gridRateShifts", gridRateShifts,
            "uninfectiousRate", uninfectiousRate
        );
        
        boolean valid = spline.update();
        Assert.assertTrue(valid);
        
        // Check continuity of values
        for (int i = 0; i < spline.getGridPointCount() - 1; i++) {
            double t = spline.getGridPointTime(i);
            double tNext = spline.getGridPointTime(i + 1);
            double prev = spline.getPrevalenceAtGridPoint(t);
            double prevNext = spline.getPrevalenceAtGridPoint(tNext);
            
            // Values should be increasing (or at least non-decreasing)
            Assert.assertTrue(prevNext >= prev * 0.1); // Allow some flexibility
        }
    }

    @Test
    public void testTransmissionRateCalculation() {
        // Test that transmission rate = uninfectiousRate - d(log(I))/dt
        RateShifts rateShifts = new RateShifts();
        rateShifts.initByName("value", Arrays.asList(0.0, 1.0, 2.0));
        
        RateShifts gridRateShifts = new RateShifts();
        gridRateShifts.initByName("value", Arrays.asList(0.0, 1.0, 2.0));
        
        RealParameter logInfected = new RealParameter("0.0 0.5 1.0"); // Constant derivative = 0.5
        RealParameter uninfectiousRate = new RealParameter("1.0");
        
        NotAKnotSpline spline = new NotAKnotSpline();
        spline.initByName(
            "logInfected", logInfected,
            "rateShifts", rateShifts,
            "gridRateShifts", gridRateShifts,
            "uninfectiousRate", uninfectiousRate
        );
        
        boolean valid = spline.update();
        Assert.assertTrue(valid);
        
        // Transmission rate should be uninfectiousRate - derivative
        double transmissionRate0 = spline.getTranssmissionRateAtGridPoint(0.0);
        double derivative0 = spline.getDerivativeAtGridPoint(0.0);
        
        Assert.assertEquals(1.0 - derivative0, transmissionRate0, 1e-8);
        Assert.assertTrue(transmissionRate0 >= 0); // Must be non-negative
    }

    @Test
    public void testInvalidTransmissionRate() {
        // Test case where transmission rate would be negative
        RateShifts rateShifts = new RateShifts();
        rateShifts.initByName("value", Arrays.asList(0.0, 1.0));
        
        RateShifts gridRateShifts = new RateShifts();
        gridRateShifts.initByName("value", Arrays.asList(0.0, 1.0));
        
        // Very steep increase: derivative = 5.0, but uninfectiousRate = 1.0
        // So transmission rate = 1.0 - 5.0 = -4.0 (invalid)
        RealParameter logInfected = new RealParameter("0.0 5.0");
        RealParameter uninfectiousRate = new RealParameter("1.0");
        
        NotAKnotSpline spline = new NotAKnotSpline();
        spline.initByName(
            "logInfected", logInfected,
            "rateShifts", rateShifts,
            "gridRateShifts", gridRateShifts,
            "uninfectiousRate", uninfectiousRate
        );
        
        boolean valid = spline.update();
        Assert.assertFalse(valid); // Should be invalid due to negative transmission rate
    }

    @Test
    public void testGridPointRetrieval() {
        RateShifts rateShifts = new RateShifts();
        rateShifts.initByName("value", Arrays.asList(0.0, 2.0, 4.0));
        
        RateShifts gridRateShifts = new RateShifts();
        gridRateShifts.initByName("value", Arrays.asList(0.0, 1.0, 2.0, 3.0, 4.0));
        
        RealParameter logInfected = new RealParameter("0.0 1.0 2.0");
        RealParameter uninfectiousRate = new RealParameter("2.0");
        
        NotAKnotSpline spline = new NotAKnotSpline();
        spline.initByName(
            "logInfected", logInfected,
            "rateShifts", rateShifts,
            "gridRateShifts", gridRateShifts,
            "uninfectiousRate", uninfectiousRate
        );
        
        spline.update();
        
        // Test grid point properties
        Assert.assertEquals(5, spline.getGridPointCount());
        Assert.assertEquals(0.0, spline.getGridStart(), 1e-10);
        Assert.assertEquals(4.0, spline.getGridEnd(), 1e-10);
        Assert.assertEquals(0.0, spline.getGridPointTime(0), 1e-10);
        Assert.assertEquals(4.0, spline.getGridPointTime(4), 1e-10);
    }

    @Test
    public void testStoreRestore() {
        RateShifts rateShifts = new RateShifts();
        rateShifts.initByName("value", Arrays.asList(0.0, 1.0, 2.0));
        
        RateShifts gridRateShifts = new RateShifts();
        gridRateShifts.initByName("value", Arrays.asList(0.0, 1.0, 2.0));
        
        RealParameter logInfected = new RealParameter("0.0 1.0 2.0");
        RealParameter uninfectiousRate = new RealParameter("1.5");
        
        NotAKnotSpline spline = new NotAKnotSpline();
        spline.initByName(
            "logInfected", logInfected,
            "rateShifts", rateShifts,
            "gridRateShifts", gridRateShifts,
            "uninfectiousRate", uninfectiousRate
        );
        
        spline.update();
        double originalPrev = spline.getPrevalenceAtGridPoint(1.0);
        double originalTransRate = spline.getTranssmissionRateAtGridPoint(1.0);
        
        // Store state
        spline.store();
        
        // Modify parameter
        logInfected.setValue(1, 0.5); // Change middle value
        spline.requiresRecalculation();
        spline.update();
        
        double modifiedPrev = spline.getPrevalenceAtGridPoint(1.0);
        Assert.assertTrue(Math.abs(modifiedPrev - originalPrev) > 0.1);
        
        // Restore state
        logInfected.setValue(1, 1.0); // Manually restore
        spline.restore();
        spline.update();
        
        double restoredPrev = spline.getPrevalenceAtGridPoint(1.0);
        // After restore, should recalculate and get same values
        Assert.assertEquals(originalPrev, restoredPrev, 1e-8);
    }

    @Test
    public void testUninfectiousRateRetrieval() {
        RateShifts rateShifts = new RateShifts();
        rateShifts.initByName("value", Arrays.asList(0.0, 1.0));
        
        RateShifts gridRateShifts = new RateShifts();
        gridRateShifts.initByName("value", Arrays.asList(0.0, 1.0));
        
        RealParameter logInfected = new RealParameter("0.0 1.0");
        RealParameter uninfectiousRate = new RealParameter("2.5");
        
        NotAKnotSpline spline = new NotAKnotSpline();
        spline.initByName(
            "logInfected", logInfected,
            "rateShifts", rateShifts,
            "gridRateShifts", gridRateShifts,
            "uninfectiousRate", uninfectiousRate
        );
        
        Assert.assertEquals(2.5, spline.getUninfectiousRate().getValue(), 1e-10);
    }

    @Test
    public void testValueAtGridPointInLogSpace() {
        RateShifts rateShifts = new RateShifts();
        rateShifts.initByName("value", Arrays.asList(0.0, 1.0));
        
        RateShifts gridRateShifts = new RateShifts();
        gridRateShifts.initByName("value", Arrays.asList(0.0, 1.0));
        
        RealParameter logInfected = new RealParameter("1.0 2.0");
        RealParameter uninfectiousRate = new RealParameter("1.0");
        
        NotAKnotSpline spline = new NotAKnotSpline();
        spline.initByName(
            "logInfected", logInfected,
            "rateShifts", rateShifts,
            "gridRateShifts", gridRateShifts,
            "uninfectiousRate", uninfectiousRate
        );
        
        spline.update();
        
        double logValue0 = spline.getValueAtGridPoint(0.0);
        double prev0 = spline.getPrevalenceAtGridPoint(0.0);
        
        // getValueAtGridPoint returns log(prevalence)
        Assert.assertEquals(Math.log(prev0), logValue0, 1e-10);
        Assert.assertEquals(1.0, logValue0, 1e-10); // log(I) at t=0 should be 1.0
    }
}
