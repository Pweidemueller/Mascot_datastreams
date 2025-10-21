package mascotdatastreams.dynamics;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.parameter.RealParameter;
import mascot.parameterdynamics.NeDynamics;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.List;

/**
 * Maps log-prevalence dynamics to an Ne(t) process using spline interpolation.
 * 
 * This class takes a NotAKnotSpline for log-prevalence interpolation and maps
 * the prevalence to coalescent effective population size Ne(t) using the transmission rate
 * formula: Ne(t) = I(t) / (c * transmission_rate(t)), where transmission_rate(t) 
 * is the derivative of log-prevalence plus the uninfectious rate.
 * 
 * Key differences from Skygrowth:
 * - Uses spline interpolation instead of piecewise-exponential
 * - Takes a pre-computed NotAKnotSpline as input rather than individual parameters
 * - Provides analytical derivatives via spline
 * - Uses precomputed grid points for efficient lookup
 */
@Description("Maps log-prevalence to Ne(t) using spline interpolation between rate shift points")
public class SplinePrevalenceToNe extends NeDynamics implements Loggable {
    
    // Required inputs
    final public Input<NotAKnotSpline> splineInput = new Input<>("spline",
            "Not-a-knot spline for log-prevalence interpolation", 
            Input.Validate.REQUIRED);
    
    
    // Optional inputs
    final public Input<RealParameter> coalescentScaleInput = new Input<>("coalescentScale",
            "Coalescent scaling constant c in Ne = I / (c * transmission_rate)", 
            Input.Validate.OPTIONAL);

    // Constants
    // TODO revisit this, since ideally we don't need to use clipping but the sampler should reject unreasonable transmissionrate values
    private static final double TSR_MIN = 1e-6;
    
    // Member variables
    private NotAKnotSpline spline;
    private RealParameter coalescentScale;
    
    boolean NesKnown = false;
    boolean returnNaN = false;
    
    
    @Override
    public void initAndValidate() {
        spline = splineInput.get();
        coalescentScale = coalescentScaleInput.get();
        
        // Set time flag for NeDynamics
        isTime = true;
        
        // Validate spline is not null
        if (spline == null) {
            throw new IllegalArgumentException("spline input is required");
        }
    }
    
    public List<String> getParameterIds() {
        return null;
    }
    
    /**
     * Gets prevalence at time t using precomputed grid points for efficiency.
     * Returns the prevalence value at the closest grid point.
     * 
     * @param t time (backward from present)
     * @return prevalence at time t
     */
    public double getPrevalenceTime(double t) {
        return spline.getPrevalenceAtGridPoint(t);
    }
    
    @Override
    public double getNeTime(double t) {
        // Get prevalence using precomputed grid points for efficiency
        double I_t = getPrevalenceTime(t);
        
        // Compute transmission rate = -dlogI/dt (forward in time)
        double transmissionRate = spline.getTranssmissionRateAtGridPoint(t);
        
        // // Clamp transmission rate to minimum value to prevent division by zero
        // transmissionRate = Math.max(transmissionRate, TSR_MIN);
        
        // Get coalescent scaling constant
        double c = coalescentScale.getArrayValue();
        
        // Compute Ne: Ne = I / (c * transmission_rate)
        return I_t / (c * transmissionRate);
    }
    
    // TODO: change to more salient recalculation criterion, e.g. just recalculate the splines and then return super.requiresRecalculation()
    @Override
    public boolean requiresRecalculation() {
        return true;  // Always recalculate when parameters change
    }
    
    // TODO do more salient isDirty return
    @Override
    public boolean isDirty() {        
        
        return true;
    }
    
    @Override
    public void store() {
        super.store();
    }

    @Override
    public void restore() {
        super.restore();
    }

    @Override
    public void init(PrintStream printStream) {
        for (int i = 0; i < spline.getGridPointCount(); i+=2) {
            printStream.print("logPrevalence_" + i + "\t");
        }
        for (int i = 0; i < spline.getGridPointCount(); i+=2) {
            printStream.print("logNe_" + i + "\t");
        }
        for (int i = 0; i < spline.getGridPointCount(); i+=20) {
            printStream.print("transmissionRate_" + i + "\t");
        }
    }

    @Override
    public void log(long l, PrintStream printStream) {
 
        for (int i = 0; i < spline.getGridPointCount(); i+=2) {
            double t = spline.getGridPointTime(i);
            double prevalence = getPrevalenceTime(t);
            printStream.print(Math.log(prevalence) + "\t");
        }
        for (int i = 0; i < spline.getGridPointCount(); i+=2) {
            double t = spline.getGridPointTime(i);
            double Ne = getNeTime(t);
            printStream.print(Math.log(Ne) + "\t");
        }
        for (int i = 0; i < spline.getGridPointCount(); i+=20) {
            double t = spline.getGridPointTime(i);
            double transmissionRate = spline.getTranssmissionRateAtGridPoint(t);
            // Clamp transmission rate to minimum value
            transmissionRate = Math.max(transmissionRate, TSR_MIN);
            printStream.print(transmissionRate + "\t");
        }
    }

    @Override
    public void close(PrintStream printStream) {
        // Nothing to close
    }
}