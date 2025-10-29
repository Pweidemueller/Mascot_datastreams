package mascotdatastreams.operators;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Operator;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

/**
 * Up-down operator that performs random walk scaling on parameters.
 * Supports both parameters already in log space and parameters in real space.
 *
 * Parameters in the "up" lists are scaled by adding a random value (log space)
 * or multiplied by a factor f (real space).
 * Parameters in the "down" lists are scaled in the opposite direction.
 *
 * This is useful for jointly scaling parameters like prevalence (log space) and
 * case count scaling factors (real space), where one should increase as the other decreases.
 */
public class UpDownLogScaleOperator extends Operator {

    final public Input<Double> scaleFactorInput = new Input<>("scaleFactor",
            "Magnitude factor used for scaling (standard deviation of random walk in log space)",
            Validate.REQUIRED);

    public Input<List<RealParameter>> upLogParametersInput = new Input<>("upLogParameter",
            "Parameters ALREADY IN LOG SPACE to scale UP (add positive random value)",
            new ArrayList<>());

    public Input<List<RealParameter>> downLogParametersInput = new Input<>("downLogParameter",
            "Parameters ALREADY IN LOG SPACE to scale DOWN (add negative random value)",
            new ArrayList<>());

    public Input<List<RealParameter>> upRealParametersInput = new Input<>("upRealParameter",
            "Parameters IN REAL SPACE to scale UP (multiply by exp(random value))",
            new ArrayList<>());

    public Input<List<RealParameter>> downRealParametersInput = new Input<>("downRealParameter",
            "Parameters IN REAL SPACE to scale DOWN (multiply by exp(-random value))",
            new ArrayList<>());

    final public Input<Boolean> optimiseInput = new Input<>("optimise",
            "Flag to indicate that the scale factor is automatically changed to achieve good acceptance rate",
            true);

    final public Input<Double> scaleUpperLimit = new Input<>("upper",
            "Upper limit of scale factor", 10.0);

    final public Input<Double> scaleLowerLimit = new Input<>("lower",
            "Lower limit of scale factor", 0.0);

    double scaleFactor;
    List<RealParameter> upLogParameters, downLogParameters;
    List<RealParameter> upRealParameters, downRealParameters;
    double upper, lower;

    @Override
    public void initAndValidate() {
        scaleFactor = scaleFactorInput.get();
        upLogParameters = upLogParametersInput.get();
        downLogParameters = downLogParametersInput.get();
        upRealParameters = upRealParametersInput.get();
        downRealParameters = downRealParametersInput.get();
        upper = scaleUpperLimit.get();
        lower = scaleLowerLimit.get();

        if (upLogParameters.isEmpty() && downLogParameters.isEmpty() &&
            upRealParameters.isEmpty() && downRealParameters.isEmpty()) {
            throw new IllegalArgumentException("At least one parameter must be specified");
        }
    }

    @Override
    public double proposal() {
        double logHR = 0.0;

        // Generate random walk step in log space using Gaussian distribution
        // scaleFactor controls the standard deviation of the walk
        double logScaler = Randomizer.nextGaussian() * scaleFactor;

        try {
            // UP parameters in LOG space: directly add logScaler
            for (RealParameter param : upLogParameters) {
                param.startEditing(this);
                for (int i = 0; i < param.getDimension(); i++) {
                    double oldValue = param.getValue(i);
                    // Directly add random walk step in log space
                    double newValue = oldValue + logScaler;
                    param.setValue(i, newValue);
                }
                // Hastings ratio contribution: logScaler for each dimension
                logHR += logScaler * param.getDimension();
            }

            // DOWN parameters in LOG space: directly subtract logScaler
            for (RealParameter param : downLogParameters) {
                param.startEditing(this);
                for (int i = 0; i < param.getDimension(); i++) {
                    double oldValue = param.getValue(i);
                    // Directly subtract random walk step in log space
                    double newValue = oldValue - logScaler;
                    param.setValue(i, newValue);
                }
                // Hastings ratio contribution: -logScaler for each dimension
                logHR -= logScaler * param.getDimension();
            }

            // UP parameters in REAL space: multiply by exp(logScaler)
            for (RealParameter param : upRealParameters) {
                param.startEditing(this);
                for (int i = 0; i < param.getDimension(); i++) {
                    double oldValue = param.getValue(i);
                    if (oldValue <= 0) {
                        throw new IllegalArgumentException(
                            "Real-space parameter " + param.getID() + " must be positive for scaling");
                    }
                    // Transform to log space, add step, transform back
                    double newValue = oldValue * Math.exp(logScaler);
                    param.setValue(i, newValue);
                }
                // Hastings ratio: logScaler for each dimension
                logHR += logScaler * param.getDimension();
            }

            // DOWN parameters in REAL space: multiply by exp(-logScaler)
            for (RealParameter param : downRealParameters) {
                param.startEditing(this);
                for (int i = 0; i < param.getDimension(); i++) {
                    double oldValue = param.getValue(i);
                    if (oldValue <= 0) {
                        throw new IllegalArgumentException(
                            "Real-space parameter " + param.getID() + " must be positive for scaling");
                    }
                    // Transform to log space, subtract step, transform back
                    double newValue = oldValue * Math.exp(-logScaler);
                    param.setValue(i, newValue);
                }
                // Hastings ratio: -logScaler for each dimension
                logHR -= logScaler * param.getDimension();
            }

        } catch (IllegalArgumentException ex) {
            return Double.NEGATIVE_INFINITY;
        }

        return logHR;
    }

    /**
     * Automatic parameter tuning based on acceptance rate
     */
    @Override
    public void optimize(final double logAlpha) {
        if (optimiseInput.get()) {
            double delta = calcDelta(logAlpha);
            delta += Math.log(scaleFactor);
            double newScaleFactor = Math.exp(delta);
            setCoercableParameterValue(newScaleFactor);
        }
    }

    @Override
    public double getCoercableParameterValue() {
        return scaleFactor;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        scaleFactor = Math.max(Math.min(value, upper), lower);
    }

    @Override
    public String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // Adjust scale factor based on acceptance rate
        final double newScaleFactor = scaleFactor * ratio;

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting scaleFactor to about " + formatter.format(newScaleFactor);
        } else if (prob > 0.40) {
            return "Try setting scaleFactor to about " + formatter.format(newScaleFactor);
        } else {
            return "";
        }
    }
}
