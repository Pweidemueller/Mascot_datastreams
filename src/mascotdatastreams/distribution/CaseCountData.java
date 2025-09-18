package mascotdatastreams.distribution;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;

@Description("Manages case count data with timestamps for each trait")
public class CaseCountData extends BEASTObject {
    final public Input<RealParameter> caseCounts = new Input<>("caseCounts", "Case counts for each observation", Validate.REQUIRED);
    final public Input<RealParameter> observationTimes = new Input<>("observationTimes", "Times of case count observations in years relative to the most recent sample", Validate.REQUIRED);
    final public Input<RealParameter> traitIndices = new Input<>("traitIndices", "Trait index for each observation (0-based)", Validate.REQUIRED);
    
    protected Double[] counts;
    protected Double[] times;
    protected int[] traits;
    protected int observationCount;
    
    @Override
    public void initAndValidate() {
        counts = caseCounts.get().getValues();
        times = observationTimes.get().getValues();
        Double[] traitValues = traitIndices.get().getValues();
        
        observationCount = counts.length;
        
        // Validate dimensions
        if (times.length != observationCount || traitValues.length != observationCount) {
            throw new IllegalArgumentException("All input arrays must have the same length");
        }
        
        // Convert trait indices to integers
        traits = new int[observationCount];
        for (int i = 0; i < observationCount; i++) {
            traits[i] = (int) Math.round(traitValues[i]);
        }
    }
    
    public int getObservationCount() {
        return observationCount;
    }
    
    public double getCaseCount(int index) {
        return counts[index];
    }
    
    public double getTime(int index) {
        return times[index];
    }
    
    public int getTrait(int index) {
        return traits[index];
    }
}
