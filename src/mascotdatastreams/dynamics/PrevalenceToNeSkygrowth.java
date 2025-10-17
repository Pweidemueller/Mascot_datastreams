package mascotdatastreams.dynamics;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import mascot.parameterdynamics.NeDynamics;
import mascot.parameterdynamics.Skygrowth;
import java.lang.reflect.Field;

/**
 * Maps log-prevalence dynamics to an Ne(t) process compatible with MASCOT callers.
 *
 * Semantics mirror mascot.parameterdynamics.Skygrowth (backward time):
 * - log I(t) is specified at breakpoints and interpolated piecewise-exponentially via a
 *   backward-time growth rate per interval.
 * - Let g_bwd be the backward-time slope of log I in the current interval.
 * - Forward-time derivative: d/dt_fwd log I = -g_bwd, so dI/dt = I * (-g_bwd).
 * - Transmission rate: transmission_rate(t) = dI/dt / I + gamma = -g_bwd + gamma.
 * - Coalescent effective population size: Ne(t) = I(t) / (c * transmission_rate(t)), with c defaulting to 2.
 */
@Description("Prevalence-to-Ne mapping with Skygrowth-style log-prevalence and RateShifts; implements NeDynamics.")
public class PrevalenceToNeSkygrowth extends NeDynamics {
    // TO DO: this should be a copy of skygrowth 
    public final Input<Skygrowth> prevalenceInput = new Input<>(
            "prevalence",
            "Mascot Skygrowth providing values that are treated as prevalence in log space (semantics: getNeTime returns exp(log I(t)))",
            Input.Validate.REQUIRED);

    public final Input<RealParameter> uninfectiousRateInput = new Input<>(
            "uninfectiousRate", "global uninfectious rate (gamma)", Input.Validate.REQUIRED);

    // Optional: coalescent scaling constant c (default 2)
    public final Input<RealParameter> coalescentScaleInput = new Input<>(
            "coalescentScale", "coalescent scaling constant c in Ne = I / (c * transmission_rate)", Input.Validate.OPTIONAL);

    private Skygrowth prevalence;
    private RealParameter uninfectiousRate;
    private RealParameter coalescentScale; // optional; if null -> use 2.0

    boolean isValid = true;

    // Numerical safety clamps
    private static final double EPS = 1e-8;
    private static final double I_MIN = 1e-12;
    private static final double I_MAX = 1e12;
    private static final double NE_MIN = 1e-6;
    private static final double NE_MAX = 1e12;

    // Step for finite-difference slope of log-prevalence (years before present)
    private static final double DT = 1e-5;

    @Override
    public void initAndValidate() {
        prevalence = prevalenceInput.get();
        uninfectiousRate = uninfectiousRateInput.get();
        coalescentScale = coalescentScaleInput.get();
        isTime = true;
    }

    @Override
    public double getNeTime(double t) {
        isValid = true;
        // Compute forward-time slope of log I(t) from Skygrowth via central difference
        double slope = getForwardSlopeAt(t);
        double I_t = prevalence.getNeTime(t);
        // Clamp I(t) to feasible numeric range
        if (I_t < I_MIN) { isValid = false; I_t = I_MIN; }
        if (I_t > I_MAX) { isValid = false; I_t = I_MAX; }
        double gamma = uninfectiousRate.getArrayValue();
        if (!(gamma >= 0.0)) {
            // Invalid gamma: clamp to nonnegative
            isValid = false;
            gamma = 0.0;
        }
        double transmission_rate = slope + gamma;
        if (!(transmission_rate > 0.0)) {
            // Negative/zero transmission rate: clamp to EPS
            isValid = false;
            transmission_rate = Math.max(transmission_rate, EPS);
        }

        double c = (coalescentScale != null) ? coalescentScale.getArrayValue() : 2.0;
        if (!(c > 0.0)) {
            // Invalid coalescent scaling constant: clamp to EPS
            isValid = false;
            c = EPS;
        }

        double Ne = I_t / (c * transmission_rate);
        if (Double.isNaN(Ne) || Double.isInfinite(Ne)) {
            isValid = false;
            Ne = (Ne >= 0.0) ? NE_MAX : NE_MIN;
        }
        if (Ne < NE_MIN) { isValid = false; Ne = NE_MIN; }
        if (Ne > NE_MAX) { isValid = false; Ne = NE_MAX; }

        return Ne;
    }

    private double getForwardSlopeAt(double t) {
        // Prefer exact behavior via reflection when possible: forward-time slope = -growth[interval]
        // and 0.0 after the last rate shift, matching Skygrowth.getNeTime() semantics.
        try {
            // access growth[]
            Field fg = prevalence.getClass().getDeclaredField("growth");
            fg.setAccessible(true);
            Object g = fg.get(prevalence);
            // access rateShifts
            Field fr = prevalence.getClass().getDeclaredField("rateShifts");
            fr.setAccessible(true);
            Object rs = fr.get(prevalence);
            if (g instanceof double[] && rs instanceof mascot.dynamics.RateShifts) {
                double[] growthArr = (double[]) g;
                mascot.dynamics.RateShifts r = (mascot.dynamics.RateShifts) rs;
                int dim = r.getDimension();
                // After the last breakpoint: Skygrowth returns a constant Ne, so slope = 0
                if (dim == 0) return 0.0;
                if (t >= r.getValue(dim - 1)) return 0.0;
                // Find interval index i such that t < r.getValue(i)
                for (int i = 0; i < dim; i++) {
                    if (t < r.getValue(i)) {
                        // forward-time slope = - backward-time growth
                        return growthArr[i];
                    }
                }
                // Fallback (should not happen due to earlier return): after last interval
                return 0.0;
            }
        } catch (Throwable ignore) { }
        return 0.0;

    }

    @Override
	public boolean requiresRecalculation() {
		// return super.requiresRecalculation();
		return true;
	}
    
    @Override
    public boolean isDirty() {
        if (uninfectiousRate != null && uninfectiousRate.isDirty(0)) return true;
        if (coalescentScale != null && coalescentScale.isDirty(0)) return true;
        return false;
    }
}
