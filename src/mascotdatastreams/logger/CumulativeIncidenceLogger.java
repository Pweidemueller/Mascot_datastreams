package mascotdatastreams.logger;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Loggable;
import beast.base.inference.parameter.RealParameter;
import mascotdatastreams.distribution.SeroprevalenceLikelihood;
import mascotdatastreams.dynamics.Spline;

import java.io.PrintStream;

/**
 * Logger for cumulative incidence values from SeroprevalenceLikelihood.
 * 
 * This logger logs cumulative incidence values at grid points of the prevalence spline.
 * It takes a SeroprevalenceLikelihood instance as input and logs the cumulative incidence
 * computed from the earliest time to each grid point.
 */
@Description("Logger for cumulative incidence values from seroprevalence likelihood")
public class CumulativeIncidenceLogger extends BEASTObject implements Loggable {
    
    public final Input<SeroprevalenceLikelihood> seroprevalenceLikelihoodInput = new Input<>(
            "seroprevalenceLikelihood", 
            "SeroprevalenceLikelihood instance to log cumulative incidence from", 
            Validate.REQUIRED);
    
    protected SeroprevalenceLikelihood seroprevalenceLikelihood;
    
    @Override
    public void initAndValidate() {
        seroprevalenceLikelihood = seroprevalenceLikelihoodInput.get();
        if (seroprevalenceLikelihood == null) {
            throw new IllegalArgumentException("CumulativeIncidenceLogger: 'seroprevalenceLikelihood' input is required");
        }
    }
    
    @Override
    public void init(PrintStream printStream) {
        // Get the prevalence spline from the seroprevalence likelihood
        Spline prevalenceSpline = seroprevalenceLikelihood.getPrevalenceSpline();
        
        // Print headers for every 10th grid point
        for (int i = 0; i < prevalenceSpline.getGridPointCount(); i += 10) {
            printStream.print("cumulativeIncidence_" + i + "\t");
        }
        for (int i = 0; i < prevalenceSpline.getGridPointCount(); i += 10) {
            printStream.print("propSeropositive_" + i + "\t");
        }

    }
    
    @Override
    public void log(long sample, PrintStream printStream) {
        // Get the prevalence spline from the seroprevalence likelihood
        Spline prevalenceSpline = seroprevalenceLikelihood.getPrevalenceSpline();
        
        // Determine earliest time (same logic as in SeroprevalenceLikelihood)
        double earliestTime = prevalenceSpline.getGridEnd();
        RealParameter earliestTimeParam = seroprevalenceLikelihood.getEarliestTimeInput();
        if (earliestTimeParam != null) {
            earliestTime = earliestTimeParam.getArrayValue();
        }
        
        // Log cumulative incidence for every 10th grid point (matching headers)
        for (int i = 0; i < prevalenceSpline.getGridPointCount(); i += 10) {
            double t = prevalenceSpline.getGridPointTime(i);
            double cumulativeIncidence = seroprevalenceLikelihood.getCumulativeIncidence(earliestTime, t);
            printStream.print(cumulativeIncidence + "\t");
        }
        for (int i = 0; i < prevalenceSpline.getGridPointCount(); i += 10) {
            double t = prevalenceSpline.getGridPointTime(i);
            double prob_seropositive = seroprevalenceLikelihood.propSeropositive(earliestTime, t);
            printStream.print(prob_seropositive + "\t");
        }
    }
    
    @Override
    public void close(PrintStream printStream) {
        // Nothing to close
    }
}

