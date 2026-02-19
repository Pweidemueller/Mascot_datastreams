package mascotdatastreams.distribution;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.parameter.RealParameter;
import mascotdatastreams.dynamics.Spline;

import java.util.List;
import java.util.Random;

/**
 * Likelihood for cumulative case count observations based on prevalence dynamics from a spline.
 * 
 * This class takes a Spline that provides log-prevalence values and computes
 * the cumulative cases by integrating the prevalence over time. The cumulative cases
 * are then used to calculate the likelihood of observed cumulative case counts.
 * 
 * The integral calculation follows the template pattern where:
 * - Prevalence I(t) is obtained from the spline
 * - Cumulative cases = integral of I(t) from start time to observation time
 * - The likelihood is computed using the cumulative cases as the mean
 */
public class SeroprevalenceLikelihood extends Distribution {
    
    // Spline input providing log-prevalence values
    public final Input<Spline> prevalenceSplineInput = new Input<>(
            "prevalenceSpline", "Spline providing log-prevalence values for interpolation", Validate.REQUIRED);

    // Observations passed directly as parameters
    public final Input<RealParameter> SeroPeopleTestedInput = new Input<>(
            "seroPeopleTested", "Observed number of people tested for seropositivity (dimension = number of observations)", Validate.REQUIRED);
    public final Input<RealParameter> SeroPeopleSeropositiveInput = new Input<>(
            "seroPeopleSeropositive", "Observed number of people seropositive (dimension = number of observations)", Validate.REQUIRED);
    public final Input<RealParameter> SeroTimesInput = new Input<>(
            "seroTimes", "Observation times corresponding 1:1 to seroPeopleSeropositive", Validate.REQUIRED);
    public final Input<RealParameter> populationSizeInput = new Input<>(
            "populationSize", "Population size", Validate.REQUIRED);
    
    final public Input<ParametricDistribution> distInput = new Input<>(
            "distribution", "Distribution used to calculate likelihood. Expected Binomial.", Validate.REQUIRED);
    
    // Optional log-space scaling factor
    public final Input<RealParameter> ScalingInput = new Input<>(
            "scaling", "Scaling factor applied to the seroprevalence proportion.", Validate.OPTIONAL);
    
    // Start time for cumulative calculation (default: time farthest away from most recent sample)
    public final Input<RealParameter> EarliestTimeInput = new Input<>(
            "EarliestTime", "Start time for cumulative case calculation (default: latest spline time)", Validate.OPTIONAL);
    
    protected RealParameter SeroPeopleTested;
    protected RealParameter SeroPeopleSeropositive;
    protected RealParameter SeroTimes;
    protected RealParameter populationSize;
    protected RealParameter Scaling;
    protected ParametricDistribution dist;
    protected Spline prevalenceSpline;

    protected boolean validated = false;
    
    @Override
    public void initAndValidate() {
        SeroPeopleTested = SeroPeopleTestedInput.get();
        SeroPeopleSeropositive = SeroPeopleSeropositiveInput.get();
        SeroTimes = SeroTimesInput.get();
        populationSize = populationSizeInput.get();
        dist = distInput.get();
        prevalenceSpline = prevalenceSplineInput.get();
        Scaling = ScalingInput.get();

        if (!(dist instanceof Binomial)) {
            throw new IllegalArgumentException("SeroprevalenceLikelihood: 'distribution' must be Binomial");
        }

        // Validate optional scaling parameter if present
        if (Scaling != null) {
            double s = Scaling.getArrayValue();
            if (!(s > 0.0)) {
                throw new IllegalArgumentException("SeroprevalenceLikelihood: 'scaling' must be > 0, got " + s);
            }
        }

        // Validate observation vectors
        if (SeroPeopleTested == null || SeroPeopleSeropositive == null || SeroTimes == null) {
            throw new IllegalArgumentException("SeroprevalenceLikelihood: 'SeroPeopleTested', 'SeroPeopleSeropositive' and 'SeroTimes' must be provided");
        }
        if (SeroPeopleTested.getDimension() != SeroPeopleSeropositive.getDimension() || SeroPeopleTested.getDimension() != SeroTimes.getDimension()) {
            throw new IllegalArgumentException("SeroprevalenceLikelihood: 'SeroPeopleTested', 'SeroPeopleSeropositive' and 'SeroTimes' must have the same dimension");
        }
        if (SeroPeopleTested.getDimension() == 0) {
            System.err.println("Warning: SeroprevalenceLikelihood received zero observations. Likelihood will be neutral (0).");
        }

        calculateLogP();
        validated = true;
    }
    
    /**
     * Computes the population-level cumulative hazard: the integral of β(t)·I(t) from EarliestTime to endTime.
     * This integrates the incidence rate assuming the entire population remains susceptible,
     * so the result equals Λ where Λ is the per-capita cumulative hazard (force of infection over time).
     * 
     * @param EarliestTime start time for integration
     * @param endTime end time for integration
     * @return cumulative hazard (integral of β·I, not a count of unique infections)
     */
    public double getPopulationCumulativeHazard(double EarliestTime, double endTime) {
        // Times are relative to the most recent sample: backward time, so if EarliestTime is before endTime, return 0
        // The integral needs to be calculated from the first infection (around tree root) to the observation time  
        if (EarliestTime <= endTime) {
            return 0.0;
        }
        
        // Grid info
        int n = prevalenceSpline.getGridPointCount();
        if (n <= 1) {
            return 0.0;
        }

        double gridMin = prevalenceSpline.getGridStart();
        double gridMax = prevalenceSpline.getGridEnd();

        // Clamp to grid range in backward time 
        double from = Math.min(EarliestTime, gridMax);
        double to = Math.max(endTime, gridMin);
        if (from <= to) {
            return 0.0;
        }

        // Find boundary indices in backward time
        int maxIdx = prevalenceSpline.getRightGridIndex(from);
        int minIdx = prevalenceSpline.getLeftGridIndex(to);
        // Integrate using trapezoids over gridpoint segments
        double sum = 0.0;

        for (int i = maxIdx; i > minIdx; i--) {
            double t1 = prevalenceSpline.getGridPointTime(i);
            double t2 = prevalenceSpline.getGridPointTime(i - 1);

            double b1 = prevalenceSpline.getTranssmissionRateAtGridPoint(t1);
            double b2 = prevalenceSpline.getTranssmissionRateAtGridPoint(t2);

            double I1 = prevalenceSpline.getPrevalenceAtGridPoint(t1);
            double I2 = prevalenceSpline.getPrevalenceAtGridPoint(t2);

            double incidence1 = b1 * I1;
            double incidence2 = b2 * I2;

            sum += (t1 - t2) * (incidence1 + incidence2) * 0.5;
        }

        // divide by population size to get per-capita cumulative hazard
        sum /= populationSize.getArrayValue();

        return sum;
    }

    /**
     * Computes the proportion of seropositive people using the current inferred scaling parameter.
     * This method automatically uses the scaling parameter from Scaling if available,
     * otherwise defaults to 1.0.
     * 
     * @param EarliestTime start time for cumulative incidence calculation
     * @param t end time for cumulative incidence calculation
     * @return proportion of seropositive people
     */
    public double propSeropositive(double EarliestTime, double t) {
        // Apply optional scaling factor (default 1.0), to adjust cumulative hazard for for more/less exposed individuals (e.g. healthcare workers vs. general population)
        // Returns the probability of an individual being seropositive at time t based on cumulative hazard and assuming Poisson process for infection with cumulative hazard as rate
        double scaling = 1.0;
        if (Scaling != null) {
            scaling = Scaling.getArrayValue();
        }
        double populationCumulativeHazard = getPopulationCumulativeHazard(EarliestTime, t);
        return 1.0 - Math.exp(-1 * scaling * populationCumulativeHazard);
    }
    
    @Override
    public double calculateLogP() {
        logP = 0.0;
        // Neutral likelihood if there are no observations
        if (SeroPeopleTested == null || SeroPeopleTested.getDimension() == 0) {
            return logP;
        }
        // Default to last grid time since everything is represented relative to most recent sample
        // we want to calculate cumulative cases from the first infection (around tree root) to the observation time  
        double EarliestTime = prevalenceSpline.getGridEnd(); 
        RealParameter EarliestTimeParam = EarliestTimeInput.get();
        if (EarliestTimeParam != null) {
            EarliestTime = EarliestTimeParam.getArrayValue();
        }

        int numObservations = SeroPeopleTested.getDimension();
        for (int i = 0; i < numObservations; i++) {
            double t = SeroTimes.getArrayValue(i);
            int n = (int) Math.round(SeroPeopleTested.getArrayValue(i));
            int x = (int) Math.round(SeroPeopleSeropositive.getArrayValue(i));
            // Probability of an individual being seropositive at time t based on cumulative incidence
            double p = propSeropositive(EarliestTime, t);
            if (Double.isNaN(p)) {
                System.err.println("Warning: Invalid 'p' value: " + p);
                return Double.NEGATIVE_INFINITY;
            }
            // Clamp to (0,1) for numerical stability and validity
            // TODO: revisit clamping 
            final double PROB_EPS = 1e-16;
            p = Math.min(1.0 - PROB_EPS, Math.max(PROB_EPS, p));

            // Validate observation bounds
            if (n < 0 || x < 0 || x > n) {
                System.err.println("Warning: Invalid observation at index " + i + ": n=" + n + ", x=" + x);
                return Double.NEGATIVE_INFINITY;
            }

            // Calculate binomial log-likelihood using the provided distribution instance
            double intervalLogP = ((Binomial) dist).logPMFForParams(x, n, p);
            if (Double.isNaN(intervalLogP)) {
                System.err.println("Warning: NaN likelihood for observation index " + i + ": n=" + n + ", x=" + x + ", p=" + p);
                return Double.NEGATIVE_INFINITY;
            }
            // if (t < 0.001) {
            //     System.out.println(">>>Start time: " + EarliestTime);
            //     System.out.println("t: " + t);
            //     System.out.println("scaling: " + scaling);
            //     System.out.println("cumulativeCases: " + cumulativeCases);
            //     System.out.println("x: " + x);
            //     System.out.println("n: " + n);
            //     System.out.println("p: " + p);
            //     System.out.println("intervalLogP: " + intervalLogP);
            // }

            logP += intervalLogP;
        }
        return logP;
    }

    public double getCumulativeIncidence(double EarliestTime, double t) {
        double populationCumulativeHazard = getPopulationCumulativeHazard(EarliestTime, t);
        return (1 - Math.exp(-1 * populationCumulativeHazard)) * populationSize.getArrayValue();
    }
    
    /**
     * Getter for prevalence spline (for use by loggers).
     * @return the prevalence spline
     */
    public Spline getPrevalenceSpline() {
        return prevalenceSpline;
    }
    
    /**
     * Getter for earliest time input (for use by loggers).
     * @return the earliest time parameter input
     */
    public RealParameter getEarliestTimeInput() {
        return EarliestTimeInput.get();
    }

    @Override
    public boolean requiresRecalculation() {
        return true;
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
        // Not implemented for likelihood classes
    }
}
