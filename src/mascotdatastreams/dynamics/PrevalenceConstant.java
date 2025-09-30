package mascotdatastreams.dynamics;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;

/**
 * Constant log-prevalence dynamics: log I(t) = c for all t.
 * Useful for tests and simple baselines without external dependencies.
 */
@Description("Constant log-prevalence dynamics: log I(t) = c for all t.")
public class PrevalenceConstant extends PrevalenceDynamics {

    public final Input<RealParameter> logPrevalenceInput = new Input<>(
            "logPrevalence", "constant log-prevalence value", Validate.REQUIRED);

    private RealParameter logPrevalence;

    @Override
    public void initAndValidate() {
        logPrevalence = logPrevalenceInput.get();
        if (logPrevalence.getDimension() != 1) {
            logPrevalence.setDimension(1);
        }
        isTime = true;
    }

    @Override
    public double getPrevalenceTime(double t) {
        return logPrevalence.getArrayValue();
    }

    @Override
    public boolean isDirty() {
        return logPrevalence.isDirty(0);
    }
}
