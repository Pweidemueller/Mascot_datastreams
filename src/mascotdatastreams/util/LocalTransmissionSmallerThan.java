package mascotdatastreams.util;


import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import mascotdatastreams.dynamics.Spline;
import mascotdatastreams.dynamics.SplinePrevalenceToNe;


@Description("returns 0 if the local (intra-deme) transmission rate is positive at every grid point, "
        + "and negative infinity if any local transmission rate is <= 0. Local transmission rate is "
        + "beta^local_i(t) = beta_i(t) - sum_{j != i} m^fw_{j->i} * I_j(t) / I_i(t), as exposed by "
        + "SplinePrevalenceToNe#getLocalTransmissionRate.")
public class LocalTransmissionSmallerThan extends Distribution {
    final public Input<SplinePrevalenceToNe> dynamicsInput = new Input<>("dynamics",
            "SplinePrevalenceToNe to check local transmission rates against", Validate.REQUIRED);


    @Override
    public void initAndValidate() {
        calculateLogP();
    }

    @Override
    public double calculateLogP() {
        SplinePrevalenceToNe dyn = dynamicsInput.get();
        Spline spline = dyn.splineInput.get();

        // Ensure the underlying spline is up to date.
        if (!spline.update()) {
            logP = Double.NEGATIVE_INFINITY;
            return Double.NEGATIVE_INFINITY;
        }

        int gridPointCount = spline.getGridPointCount();
        for (int i = 0; i < gridPointCount; i++) {
            double t = spline.getGridPointTime(i);
            if (dyn.getLocalTransmissionRate(t) <= 0) {
                logP = Double.NEGATIVE_INFINITY;
                return Double.NEGATIVE_INFINITY;
            }
        }

        logP = 0.0;
        return 0.0;
    }


    @Override
    public void sample(State state, Random random) {
        // do nothing
    }

    @Override
    public List<String> getConditions() {
        return new ArrayList<>();
    }

    @Override
    public List<String> getArguments() {
        return new ArrayList<>();
    }
}
