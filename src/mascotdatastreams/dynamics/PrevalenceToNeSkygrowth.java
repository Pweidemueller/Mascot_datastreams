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
 * - Transmission rate: transmission_rate(t) = dI/dt / I + gamma = -g_bwd + gamma.
 * - Coalescent effective population size: Ne(t) = I(t) / (c * transmission_rate(t)), with c defaulting to 2.
 */
@Description("Prevalence-to-Ne mapping with Skygrowth-style log-prevalence and RateShifts; implements NeDynamics.")
public class PrevalenceToNeSkygrowth extends NeDynamics {

    public final Input<RealParameter> logPrevalenceInput = new Input<>(
            "logPrevalence", "log-prevalence values at breakpoints", Input.Validate.REQUIRED);

    public final Input<RateShifts> rateShiftsInput = new Input<>(
            "rateShifts",
            "time breakpoints provided via RateShifts; values are specified as fractions of the tree root height (0..1) and resolved to absolute times (ascending, years before most recent sample)",
            Input.Validate.REQUIRED);

    public final Input<RealParameter> uninfectiousRateInput = new Input<>(
            "uninfectiousRate", "global uninfectious rate (gamma)", Input.Validate.REQUIRED);

    // Optional: coalescent scaling constant c (default 2)
    public final Input<RealParameter> coalescentScaleInput = new Input<>(
            "coalescentScale", "coalescent scaling constant c in Ne = I / (c * transmission_rate)", Input.Validate.OPTIONAL);

    private RealParameter logPrevalence;
    private RateShifts rateShifts;
    private RealParameter uninfectiousRate;
    private RealParameter coalescentScale; // optional; if null -> use 2.0

    private double[] growth;        // forward-time slopes of log I per interval
    private double[] growth_stored;

    // Precomputed Ne at control points (in log-space) and their interval growth
    private double[] logNeCtrl;            // length = rateShifts.getDimension() + 1
    private double[] logNeCtrl_stored;
    private double[] neGrowth;             // length = rateShifts.getDimension()
    private double[] neGrowth_stored;

    boolean isValid = true;

    // Numerical safety clamps
    private static final double EPS = 1e-8;
    private static final double I_MIN = 1e-12;
    private static final double I_MAX = 1e12;
    private static final double NE_MIN = 1e-6;
    private static final double NE_MAX = 1e12;

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
        recalcNeBreakpoints();
        isTime = true;
    }

    @Override
    public double getNeTime(double t) {
        // Interpolate on Ne level (log-space), akin to mascot.parameterdynamics.Skygrowth
        int interval = getIntervalNr(t);
        double Ne;
        if (interval >= rateShifts.getDimension()) {
            Ne = Math.exp(logNeCtrl[logNeCtrl.length - 1]);
        } else {
            double timediff = t;
            if (interval > 0) timediff -= rateShifts.getValue(interval - 1);
            double logNe_t = logNeCtrl[interval] - neGrowth[interval] * timediff;
            Ne = Math.exp(logNe_t);
            
        }
        if (Ne < NE_MIN) Ne = NE_MIN;
        if (Ne > NE_MAX) Ne = NE_MAX;
        return Ne;
    }

    private int getIntervalNr(double t) {
        for (int i = 0; i < rateShifts.getDimension(); i++)
            if (t < rateShifts.getValue(i))
                return i;
        return rateShifts.getDimension();
    }

    private void recalcGrowth() {
        // Growth is calculated forwards in time (t is expressed backwards in time)
        growth = new double[rateShifts.getDimension()];
        double curr_time = 0.0;
        for (int i = 1; i < logPrevalence.getDimension(); i++) {
            double dt = rateShifts.getValue(i - 1) - curr_time;
            growth[i - 1] = (logPrevalence.getArrayValue(i - 1) - logPrevalence.getArrayValue(i)) / dt;
            curr_time = rateShifts.getValue(i - 1);
        }
    }

    // Precompute log Ne at control points using right-hand interval slope, then interval growth on log Ne
    private void recalcNeBreakpoints() {
        int nIntervals = rateShifts.getDimension();
        int nCtrl = nIntervals + 1;
        if (logNeCtrl == null || logNeCtrl.length != nCtrl) {
            logNeCtrl = new double[nCtrl];
        }
        if (neGrowth == null || neGrowth.length != nIntervals) {
            neGrowth = new double[nIntervals];
        }

        isValid = true;
        double gamma = uninfectiousRate.getArrayValue();
        if (!(gamma >= 0.0)) { isValid = false; gamma = 0.0; }
        double c = (coalescentScale != null) ? coalescentScale.getArrayValue() : 2.0;
        if (!(c > 0.0)) { isValid = false; c = EPS; }

        // Compute log Ne values at control points
        for (int i = 0; i < nCtrl; i++) {
            double logI = logPrevalence.getArrayValue(i);
            // clamp I in log-space
            if (logI < Math.log(I_MIN)) { isValid = false; logI = Math.log(I_MIN); }
            if (logI > Math.log(I_MAX)) { isValid = false; logI = Math.log(I_MAX); }

            double g_point = (i < nIntervals) ? growth[i] : growth[nIntervals - 1];
            double transmission = g_point + gamma;
            if (!(transmission > 0.0)) { isValid = false; transmission = Math.max(transmission, EPS); }

            double logNe = logI - Math.log(c) - Math.log(transmission);
            // clamp Ne bounds in log-space
            if (logNe < Math.log(NE_MIN)) { isValid = false; logNe = Math.log(NE_MIN); }
            if (logNe > Math.log(NE_MAX)) { isValid = false; logNe = Math.log(NE_MAX); }
            logNeCtrl[i] = logNe;
        }

        // Compute growth on log Ne per interval
        double curr_time = 0.0;
        for (int i = 1; i < nCtrl; i++) {
            double dt = rateShifts.getValue(i - 1) - curr_time;
            neGrowth[i - 1] = (logNeCtrl[i - 1] - logNeCtrl[i]) / dt;
            curr_time = rateShifts.getValue(i - 1);
        }
    }

    @Override
    public boolean requiresRecalculation() {
        recalcGrowth();
        recalcNeBreakpoints();
        return super.requiresRecalculation();
    }

    @Override
    public void store() {
        growth_stored = new double[growth.length];
        System.arraycopy(growth, 0, growth_stored, 0, growth.length);
        if (logNeCtrl != null) {
            logNeCtrl_stored = new double[logNeCtrl.length];
            System.arraycopy(logNeCtrl, 0, logNeCtrl_stored, 0, logNeCtrl.length);
        }
        if (neGrowth != null) {
            neGrowth_stored = new double[neGrowth.length];
            System.arraycopy(neGrowth, 0, neGrowth_stored, 0, neGrowth.length);
        }
        super.store();
    }

    @Override
    public void restore() {
        System.arraycopy(growth_stored, 0, growth, 0, growth_stored.length);
        if (logNeCtrl_stored != null) {
            logNeCtrl = new double[logNeCtrl_stored.length];
            System.arraycopy(logNeCtrl_stored, 0, logNeCtrl, 0, logNeCtrl_stored.length);
        }
        if (neGrowth_stored != null) {
            neGrowth = new double[neGrowth_stored.length];
            System.arraycopy(neGrowth_stored, 0, neGrowth, 0, neGrowth_stored.length);
        }
        super.restore();
    }

    @Override
    public void recalculate() {
        recalcGrowth();
        recalcNeBreakpoints();
    }

    @Override
    public boolean isDirty() {
        if (logPrevalence.isDirty(0)) return true;
        if (uninfectiousRate != null && uninfectiousRate.isDirty(0)) return true;
        if (coalescentScale != null && coalescentScale.isDirty(0)) return true;
        return false;
    }
}
