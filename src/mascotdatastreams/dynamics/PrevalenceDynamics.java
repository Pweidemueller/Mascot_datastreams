package mascotdatastreams.dynamics;

import beast.base.inference.CalculationNode;

/**
 * Abstract base for prevalence dynamics per deme.
 * Implementations should provide log-prevalence trajectories over time.
 */
public abstract class PrevalenceDynamics extends CalculationNode {

    public boolean isTime;

    @Override
    public void initAndValidate() {
    }

    /**
     * Recalculate internal state when parameters change.
     */
    public void recalculate() {}

    /**
     * Return the log-prevalence log I(t) at time t (years before most recent sample).
     */
    public double getPrevalenceTime(double t) {
        throw new IllegalArgumentException("Function not implemented. Class of prevalence function not correctly recognized");
    }

    /**
     * Optional: prevalence by interval index.
     */
    public double getPrevalenceInterval(int i) {
        throw new IllegalArgumentException("Function not implemented. Class of prevalence function not correctly recognized");
    }

    public void setNrIntervals(int intervals) {}

    public boolean isDirty() {
        return true;
    }
}
