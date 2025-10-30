package mascotdatastreams.util;


import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.*;
import beast.base.inference.parameter.RealParameter;
import mascotdatastreams.dynamics.Spline;


@Description("returns 0 if all transmission rates are positive and negative infinity if any transmission rate is <= 0")
public class TransmissionSmallerThan extends Distribution {
    final public Input<Spline> splineInput = new Input<>("spline", "Spline object to check transmission rates", Validate.REQUIRED);


    @Override
    public void initAndValidate() {
        calculateLogP();
    }

    @Override
    public double calculateLogP() {
        Spline spline = splineInput.get();
        
        // Update the spline to ensure transmission rates are calculated
        if (!spline.update()) {
            logP = Double.NEGATIVE_INFINITY;
            return Double.NEGATIVE_INFINITY;
        }
        
        // Check all transmission rates at grid points
        int gridPointCount = spline.getGridPointCount();
        for (int i = 0; i < gridPointCount; i++) {
            double transmissionRate = spline.getTranssmissionRateAtGridPoint(spline.getGridPointTime(i));
            if (transmissionRate <= 0) {
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
        List<String> conditions = new ArrayList<>();
        return conditions;
    }

    @Override
    public List<String> getArguments() {
        List<String> arguments = new ArrayList<>();

        return arguments;
    }
}
