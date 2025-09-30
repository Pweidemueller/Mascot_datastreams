package mascotdatastreams.dynamics;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import mascot.dynamics.RateShifts;
import mascot.parameterdynamics.NeDynamics;

/**
 * Maps log-prevalence dynamics to an Ne(t) process compatible with MASCOT callers.
 *
 * Semantics mirror mascot.parameterdynamics.Skygrowth (backward time):
 * - log I(t) is specified at breakpoints and interpolated piecewise-exponentially via a
 *   backward-time growth rate per interval.
 * - Let g_bwd be the backward-time slope of log I in the current interval.
 * - Forward-time derivative: d/dt_fwd log I = -g_bwd, so dI/dt = I * (-g_bwd).
 * - Transmission rate: beta(t) = dI/dt / I + gamma = -g_bwd + gamma.
 * - Coalescent effective population size: Ne(t) = I(t) / (c * beta(t)), with c defaulting to 2.
 */
@Description("Prevalence-to-Ne mapping with Skygrowth-style log-prevalence and RateShifts; implements NeDynamics.")
public class PrevalenceToNeSkygrowth extends NeDynamics {

    public final Input<RealParameter> logPrevalenceInput = new Input<>(
            "logPrevalence", "log-prevalence values at breakpoints", Input.Validate.REQUIRED);

    public final Input<RateShifts> rateShiftsInput = new Input<>(
            "rateShifts", "time breakpoints (ascending, years before most recent sample)", Input.Validate.REQUIRED);

    public final Input<RealParameter> uninfectiousRateInput = new Input<>(
            "uninfectiousRate", "global uninfectious rate (gamma)", Input.Validate.REQUIRED);

    // Optional: coalescent scaling constant c (default 2)
    public final Input<RealParameter> coalescentScaleInput = new Input<>(
            "coalescentScale", "coalescent scaling constant c in Ne = I / (c * beta)", Input.Validate.OPTIONAL);

    private RealParameter logPrevalence;
    private RateShifts rateShifts;
    private RealParameter uninfectiousRate;
    private RealParameter coalescentScale; // optional; if null -> use 2.0

    private double[] growth;        // backward-time slopes of log I per interval
    private double[] growth_stored;

    @Override
    public void initAndValidate() {
        logPrevalence = logPrevalenceInput.get();
        rateShifts = rateShiftsInput.get();
        uninfectiousRate = uninfectiousRateInput.get();
        coalescentScale = coalescentScaleInput.get();

        // Ensure dimension matches number of intervals + 1
        logPrevalence.setDimension(rateShifts.getDimension() + 1);
        growth = new double[rateShifts.getDimension()];
        recalcGrowth();
        isTime = true;
    }

    @Override
    public double getNeTime(double t) {
        int interval = getIntervalNr(t);

        // Compute log I(t)
        double logI_t;
        if (interval >= rateShifts.getDimension()) {
            logI_t = logPrevalence.getArrayValue(logPrevalence.getDimension() - 1);
        } else {
            double timediff = t;
            if (interval > 0) timediff -= rateShifts.getValue(interval - 1);
            logI_t = logPrevalence.getArrayValue(interval) - growth[interval] * timediff;
        }
        double I_t = Math.exp(logI_t);
        if (!(I_t > 0.0)) {
            throw new IllegalArgumentException("PrevalenceToNeSkygrowth: Non-positive prevalence encountered at t=" + t + ".");
        }

        // Backward-time slope in the interval (or 0 after last interval)
        double g_bwd = (interval >= rateShifts.getDimension()) ? 0.0 : growth[interval];

        // Forward-time beta = -g_bwd + gamma
        double gamma = uninfectiousRate.getArrayValue();
        if (!(gamma >= 0.0)) {
            throw new IllegalArgumentException("PrevalenceToNeSkygrowth: uninfectiousRate (gamma) must be >= 0, got " + gamma);
        }
        double beta = -g_bwd + gamma;
        if (!(beta > 0.0)) {
            throw new IllegalArgumentException("PrevalenceToNeSkygrowth: implied transmission rate beta <= 0 at t=" + t +
                    ". Ensure log-prevalence smoothing and gamma produce positive beta.");
        }

        double c = (coalescentScale != null) ? coalescentScale.getArrayValue() : 2.0;
        if (!(c > 0.0)) {
            throw new IllegalArgumentException("PrevalenceToNeSkygrowth: coalescentScale c must be > 0, got " + c);
        }

        return I_t / (c * beta);
    }

    private int getIntervalNr(double t) {
        for (int i = 0; i < rateShifts.getDimension(); i++)
            if (t < rateShifts.getValue(i))
                return i;
        return rateShifts.getDimension();
    }

    private void recalcGrowth() {
        growth = new double[rateShifts.getDimension()];
        double curr_time = 0.0;
        for (int i = 1; i < logPrevalence.getDimension(); i++) {
            double dt = rateShifts.getValue(i - 1) - curr_time;
            growth[i - 1] = (logPrevalence.getArrayValue(i - 1) - logPrevalence.getArrayValue(i)) / dt;
            curr_time = rateShifts.getValue(i - 1);
        }
    }

    @Override
    public boolean requiresRecalculation() {
        recalcGrowth();
        return super.requiresRecalculation();
    }

    @Override
    public void store() {
        growth_stored = new double[growth.length];
        System.arraycopy(growth, 0, growth_stored, 0, growth.length);
        super.store();
    }

    @Override
    public void restore() {
        System.arraycopy(growth_stored, 0, growth, 0, growth_stored.length);
        super.restore();
    }

    @Override
    public void recalculate() {
        recalcGrowth();
    }

    @Override
    public boolean isDirty() {
        if (logPrevalence.isDirty(0)) return true;
        if (uninfectiousRate != null && uninfectiousRate.isDirty(0)) return true;
        if (coalescentScale != null && coalescentScale.isDirty(0)) return true;
        return false;
    }
}
