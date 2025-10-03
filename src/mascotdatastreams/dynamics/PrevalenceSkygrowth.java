package mascotdatastreams.dynamics;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import mascot.dynamics.RateShifts;

/**
 * Piecewise-exponential interpolation of log-prevalence over time, mirroring Mascot's Skygrowth for Ne.
 * The control points are in log-space and defined at rate-shift breakpoints. Within each interval,
 * log-prevalence changes linearly with slope equal to the backward-time growth rate.
 */
@Description("Prevalence function with log-scale values at time breakpoints, interpolated piecewise exponentially.")
public class PrevalenceSkygrowth extends PrevalenceDynamics {

    public final Input<RealParameter> logPrevalenceInput = new Input<>(
            "logPrevalence", "log-prevalence values at breakpoints", Input.Validate.REQUIRED);
    public final Input<RateShifts> rateShiftsInput = new Input<>(
            "rateShifts",
            "time breakpoints provided via RateShifts; values are specified as fractions of the tree root height (0..1) and resolved to absolute times (ascending, years before most recent sample)",
            Input.Validate.REQUIRED);

    RealParameter logPrevalence;
    RateShifts rateShifts;

    boolean known = false;
    double[] growth;          // per-interval backward-time growth rate of log I
    double[] growth_stored;

    @Override
    public void initAndValidate() {
        logPrevalence = logPrevalenceInput.get();
        rateShifts = rateShiftsInput.get();
        logPrevalence.setDimension(rateShifts.getDimension() + 1);
        growth = new double[rateShifts.getDimension()];
        recalc();
        isTime = true;
    }

    @Override
    public double getPrevalenceTime(double t) {
        int interval = getIntervalNr(t);
        if (interval >= rateShifts.getDimension()) {
            return logPrevalence.getArrayValue(logPrevalence.getDimension() - 1);
        }
        double timediff = t;
        if (interval > 0)
            timediff -= rateShifts.getValue(interval - 1);
        return logPrevalence.getArrayValue(interval) - growth[interval] * timediff;
    }

    private int getIntervalNr(double t) {
        for (int i = 0; i < rateShifts.getDimension(); i++)
            if (t < rateShifts.getValue(i))
                return i;
        return rateShifts.getDimension();
    }

    private void recalc() {
        growth = new double[rateShifts.getDimension()];
        double curr_time = 0.0;
        for (int i = 1; i < logPrevalence.getDimension(); i++) {
            growth[i - 1] = (logPrevalence.getArrayValue(i - 1) - logPrevalence.getArrayValue(i)) /
                    (rateShifts.getValue(i - 1) - curr_time);
            curr_time = rateShifts.getValue(i - 1);
        }
        known = true;
    }

    @Override
    public boolean requiresRecalculation() {
        recalc();
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
        // Delegate to recalc when BEAST requests
        recalc();
    }

    @Override
    public boolean isDirty() {
        if (logPrevalence.isDirty(0))
            return true;
        return false;
    }
}
