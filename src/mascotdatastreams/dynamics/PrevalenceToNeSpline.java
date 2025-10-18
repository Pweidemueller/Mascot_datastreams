package mascotdatastreams.dynamics;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import mascot.parameterdynamics.NeDynamics;
import mascot.dynamics.RateShifts;

/**
 * Maps log-prevalence dynamics to an Ne(t) process using spline interpolation.
 * 
 * This class takes log-prevalence values at rate shift points and uses cubic spline
 * interpolation to create smooth prevalence trajectories. The prevalence is then
 * mapped to coalescent effective population size Ne(t) using the transmission rate
 * formula: Ne(t) = I(t) / (c * transmission_rate(t)), where transmission_rate(t) 
 * is the derivative of log-prevalence plus the uninfectious rate.
 * 
 * Key differences from Skygrowth:
 * - Uses spline interpolation instead of piecewise-exponential
 * - Stores log-prevalence directly rather than growth rates
 * - Provides analytical derivatives via spline
 * - Configurable number of grid points for interpolation efficiency
 */
@Description("Maps log-prevalence to Ne(t) using spline interpolation between rate shift points")
public class PrevalenceToNeSpline extends NeDynamics {
    
    // Required inputs
    final public Input<RealParameter> logInfectedInput = new Input<>("logInfected",
            "Log-prevalence values at rate shift points (dimension = rateShifts + 1)", 
            Input.Validate.REQUIRED);
    
    final public Input<RateShifts> rateShiftsInput = new Input<>("rateShifts",
            "Time points for rate shifts (backward time from present)", 
            Input.Validate.REQUIRED);
    
    final public Input<RealParameter> uninfectiousRateInput = new Input<>("uninfectiousRate",
            "Rate at which individuals become uninfectious (gamma)", 
            Input.Validate.REQUIRED);
    
    // Optional inputs
    final public Input<RealParameter> coalescentScaleInput = new Input<>("coalescentScale",
            "Coalescent scaling constant c in Ne = I / (c * transmission_rate)", 
            Input.Validate.OPTIONAL);
    
    final public Input<Integer> numGridPointsInput = new Input<>("numGridPoints",
            "Number of grid points for spline interpolation (default: 1000)", 
            Input.Validate.OPTIONAL);

            
    
    // Member variables
    private RealParameter logInfected;
    private RateShifts rateShifts;
    private RealParameter uninfectiousRate;
    private RealParameter coalescentScale;
    private Integer numGridPoints;
    
    private NotAKnotSpline spline;
    private boolean splineValid = false;
    
    // TODO remove clamping where possible since ideally this should be handled by the sampler.
    // Numerical safety bounds
    private static final double EPS = 1e-8;
    private static final double I_MIN = 1e-12;
    private static final double I_MAX = 1e12;
    private static final double NE_MIN = 1e-6;
    private static final double NE_MAX = 1e12;
    
    // Default coalescent scaling constant
    private static final double DEFAULT_COALESCENT_SCALE = 2.0;
    
    @Override
    public void initAndValidate() {
        logInfected = logInfectedInput.get();
        rateShifts = rateShiftsInput.get();
        uninfectiousRate = uninfectiousRateInput.get();
        coalescentScale = coalescentScaleInput.get();
        numGridPoints = numGridPointsInput.get();
        
        // Set time flag for NeDynamics
        isTime = true;
        
        // Validate dimensions
        int nShifts = rateShifts.getDimension();
        int nLogValues = logInfected.getDimension();
        
        if (nLogValues != nShifts + 1) {
            throw new IllegalArgumentException("logInfected dimension (" + nLogValues + 
                ") must equal rateShifts dimension (" + nShifts + ") + 1");
        }

        // Validate uninfectious rate is positive
        double gamma = uninfectiousRate.getArrayValue();
        if (gamma < 0.0) {
            throw new IllegalArgumentException("uninfectiousRate must be non-negative");
        }
        
        // Validate coalescent scale if provided
        if (coalescentScale != null) {
            double c = coalescentScale.getArrayValue();
            if (c <= 0.0) {
                throw new IllegalArgumentException("coalescentScale must be positive");
            }
        }
        
        // Validate numGridPoints if provided
        if (numGridPoints != null) {
            if (numGridPoints < 2) {
                throw new IllegalArgumentException("numGridPoints must be >= 2");
            }
        }
        
        // Calculate spline coefficients and grid points
        updateSpline();
    }
    
    /**
     * Updates the spline coefficients from current input values.
     */
    private void updateSpline() {
        int nShifts = rateShifts.getDimension();
        int nPoints = nShifts + 1;
        
        double[] times = new double[nPoints];
        double[] values = new double[nPoints];
        
        // Build time points: rate shifts + present time (0.0)
        times[0] = 0.0;  // Present time
        values[0] = logInfected.getArrayValue(0);  // Log-prevalence at present
        
        for (int i = 0; i < nShifts; i++) {
            times[i + 1] = rateShifts.getValue(i);
            System.out.println("times[" + (i + 1) + "]: " + times[i + 1]);
            values[i + 1] = logInfected.getArrayValue(i + 1);
            System.out.println("values[" + (i + 1) + "]: " + values[i + 1]);
        }
        
        // Create spline with specified number of grid points (default: 1000)
        int gridPoints = (numGridPoints != null) ? numGridPoints : 1000;
        spline = new NotAKnotSpline(times, values, gridPoints);
        splineValid = true;
    }
    
    /**
     * Gets prevalence at time t using precomputed grid points for efficiency.
     * Returns the prevalence value at the closest grid point.
     * 
     * @param t time (backward from present)
     * @return prevalence at time t
     */
    public double getPrevalenceTime(double t) {
        if (!splineValid) {
            updateSpline();
        }
        double logI = spline.getValueAtGridPoint(t);
        double I = Math.exp(logI);
        return I;
    }

    public double getPrevalenceTimeExact(double t) {
        if (!splineValid) {
            updateSpline();
        }
        double logI = spline.getValue(t);
        double I = Math.exp(logI);
        return I;
    }
    
    @Override
    public double getNeTime(double t) {
        // Get derivative of log-prevalence using precomputed grid points for consistency
        if (!splineValid) {
            updateSpline();
        }
        // Get prevalence using precomputed grid points for efficiency
        double I_t = getPrevalenceTime(t);
        
        double dlogI_dt = spline.getDerivativeAtGridPoint(t);
        
        // Get uninfectious rate
        double gamma = uninfectiousRate.getArrayValue();
        System.out.println("gamma: " + gamma);
        System.out.println("dlogI_dt: " + dlogI_dt);
        
        // Compute transmission rate = dI/dt / I + gamma
        // logI = spline polynomial -> dlogI/dt = spline derivative
        // I = exp(logI) -> dI/dt = I * dlogI/dt
        // so transmission rate = dI/dt / I + gamma = (I * dlogI/dt) / I + gamma = dlogI/dt + gamma
        // t is backward in time, so dlogI/dt is the backward-time slope of log-prevalence.
        // so the forward in time slope is -dlogI/dt.
        // so the forward in time transmission rate is -dlogI/dt + gamma.
        double transmissionRate = -1 * dlogI_dt + gamma;
        // TODO revisit this clamping since it might interfere with inference/convergence.
        if (transmissionRate <= 0.0) {
            transmissionRate = EPS;  // Safety clamp
        }
        System.out.println("transmissionRate: " + transmissionRate);
        // Get coalescent scaling constant
        double c = (coalescentScale != null) ? coalescentScale.getArrayValue() : DEFAULT_COALESCENT_SCALE;
        
        // Compute Ne: Ne = I / (c * transmission_rate)
        double Ne = I_t / (c * transmissionRate);
        System.out.println("Ne: " + Ne);
        // Safety checks and clamping
        if (Double.isNaN(Ne) || Double.isInfinite(Ne)) {
            Ne = (Ne >= 0.0) ? NE_MAX : NE_MIN;
        }
        if (Ne < NE_MIN) Ne = NE_MIN;
        if (Ne > NE_MAX) Ne = NE_MAX;
        
        return Ne;
    }
    
    // TODO: change to more salient recalculation criterion, e.g. just recalculate the splines and then return super.requiresRecalculation()
    @Override
    public boolean requiresRecalculation() {
        return true;  // Always recalculate when parameters change
    }
    
    @Override
    public boolean isDirty() {
        // Check if any input parameters are dirty
        if (logInfected.isDirty(0)) return true;
        if (uninfectiousRate.isDirty(0)) return true;
        if (coalescentScale != null && coalescentScale.isDirty(0)) return true;
        
        return false;
    }
}