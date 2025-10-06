package mascotdatastreams.dynamics;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import mascot.parameterdynamics.NeDynamics;

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

    public final Input<PrevalenceSkygrowth> prevalenceInput = new Input<>(
            "prevalence",
            "PrevalenceSkygrowth providing log-prevalence interpolation and interval slopes",
            Input.Validate.REQUIRED);

    public final Input<RealParameter> uninfectiousRateInput = new Input<>(
            "uninfectiousRate", "global uninfectious rate (gamma)", Input.Validate.REQUIRED);

    // Optional: coalescent scaling constant c (default 2)
    public final Input<RealParameter> coalescentScaleInput = new Input<>(
            "coalescentScale", "coalescent scaling constant c in Ne = I / (c * transmission_rate)", Input.Validate.OPTIONAL);

    private PrevalenceSkygrowth prevalence;
    private RealParameter uninfectiousRate;
    private RealParameter coalescentScale; // optional; if null -> use 2.0

    boolean isValid = true;

    // Numerical safety clamps
    private static final double EPS = 1e-8;
    private static final double I_MIN = 1e-12;
    private static final double I_MAX = 1e12;
    private static final double NE_MIN = 1e-6;
    private static final double NE_MAX = 1e12;

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
        // Compute log I(t) and forward-time slope using the provided prevalence dynamics
        double logI_t = prevalence.getPrevalenceTime(t);
        double g_fwd = prevalence.getForwardSlopeAt(t);
        double I_t = Math.exp(logI_t);
        // Infeasible I(t): clamp to small positive
        if (I_t < I_MIN) { isValid = false; I_t = I_MIN; }
        if (I_t > I_MAX) { isValid = false; I_t = I_MAX; }

        // Forward-time transmission_rate = g_fwd + gamma
        double gamma = uninfectiousRate.getArrayValue();
        if (!(gamma >= 0.0)) {
            // Invalid gamma: clamp to nonnegative
            isValid = false;
            gamma = 0.0;
        }
        double transmission_rate = g_fwd + gamma;
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

    // Interval bookkeeping is delegated to the prevalence dynamics.

    @Override
    public boolean isDirty() {
        if (prevalence != null && prevalence.isDirty()) return true;
        if (uninfectiousRate != null && uninfectiousRate.isDirty(0)) return true;
        if (coalescentScale != null && coalescentScale.isDirty(0)) return true;
        return false;
    }
}
