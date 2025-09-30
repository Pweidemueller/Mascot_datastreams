package mascotdatastreams.dynamics;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * Container for per-deme prevalence dynamics (log-prevalence trajectories).
 */
public class PrevalenceDynamicsList extends CalculationNode {
    public Input<List<PrevalenceDynamics>> prevalenceDynamicsInput = new Input<>(
            "prevalence", "log-prevalence dynamics per deme", new ArrayList<>(), Validate.REQUIRED);

    List<PrevalenceDynamics> prevalenceDynamics;

    // Optional mapping when built via BEAUti
    public HashMap<String, Integer> traitToType = new HashMap<>();

    public int nrIntervals = 1;

    public PrevalenceDynamicsList() {
        prevalenceDynamics = prevalenceDynamicsInput.get();
    }

    @Override
    public void initAndValidate() {
        prevalenceDynamics = prevalenceDynamicsInput.get();
    }

    public int size() {
        return prevalenceDynamicsInput.get().size();
    }

    public PrevalenceDynamics get(int index) {
        if (prevalenceDynamicsInput.get() == null)
            return null;
        return prevalenceDynamicsInput.get().get(index);
    }

    public void add(PrevalenceDynamics dyn) {
        prevalenceDynamicsInput.get().add(dyn);
    }

    public void clear() {
        prevalenceDynamicsInput.get().clear();
    }
}
