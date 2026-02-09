package mascotdatastreams.logger;

import beast.base.core.Citation;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.IntervalType;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.BooleanParameter;
import mascot.distribution.Mascot;
import mascot.distribution.StructuredTreeIntervals;
import mascot.dynamics.Dynamics;
import mascot.ode.Euler2ndOrderTransitions;
import mascot.ode.MascotODEUpDown;
import mascotdatastreams.distribution.MascotLogPflag;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;

import java.io.PrintStream;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;


/**
 * Structured tree logger; accepts either mascot.distribution.Mascot or mascotdatastreams.distribution.MascotLogPflag.
 *
 * @author Nicola Felix Mueller (nicola.felix.mueller@gmail.com)
 */
@Citation(	"Nicola F. Müller, David A. Rasmussen, Tanja Stadler (2018)\n"+
			"  MASCOT: parameter and state inference under the marginal\n"+
			"  structured coalescent approximation\n"+
			"  Bioinformatics, , bty406, https://doi.org/10.1093/bioinformatics/bty406")
public class StructuredTreeLogger extends Tree implements Loggable {

	public Input<Mascot> mascotInput = new Input<>("mascot", "Mascot tree distribution (use this or mascotLogPflag)", Input.Validate.OPTIONAL);
	public Input<MascotLogPflag> mascotLogPflagInput = new Input<>("mascotLogPflag", "MascotLogPflag tree distribution (use this or mascot)", Input.Validate.OPTIONAL);

	public Input<Double> epsilonInput = new Input<>("epsilon", "step size for the RK4 integration",0.00001);
	public Input<Double> maxStepInput = new Input<>("maxStep", "step size for the RK4 integration", Double.POSITIVE_INFINITY);
	public Input<Double> stepSizeInput = new Input<>("stepSize", "step size for the RK4 integration");

	public Input<Boolean> useUpDown = new Input<>("upDown", "if up down algorithm is to use for the node state calculation", true);


    public Input<BranchRateModel.Base> clockModelInput = new Input<BranchRateModel.Base>("branchratemodel", "rate to be logged with branches of the tree");
    public Input<List<Function>> parameterInput = new Input<List<Function>>("metadata", "meta data to be logged with the tree nodes", new ArrayList<>());
    public Input<Boolean> maxStateInput = new Input<Boolean>("maxState", "report branch lengths as substitutions (branch length times clock rate for the branch)", false);
    public Input<BooleanParameter> conditionalStateProbsInput = new Input<BooleanParameter>("conditionalStateProbs", "report branch lengths as substitutions (branch length times clock rate for the branch)");
    public Input<Boolean> substitutionsInput = new Input<Boolean>("substitutions", "report branch lengths as substitutions (branch length times clock rate for the branch)", false);
    public Input<Integer> decimalPlacesInput = new Input<Integer>("dp", "the number of decimal places to use writing branch lengths and rates, use -1 for full precision (default = full precision)", -1);


    protected boolean someMetaDataNeedsLogging;
    protected boolean substitutions = false;
    protected boolean takeMax = true;
    protected boolean conditionals = true;

    protected boolean updown = true;

    protected DecimalFormat df;
    protected String type;

    protected int states;
    protected boolean[] used;
    protected boolean report;
    /** Either Mascot or MascotLogPflag. */
    protected Object mascotSource;

	TreeInterface tree;


	private StructuredTreeIntervals getStructuredTreeIntervals() {
		if (mascotSource instanceof Mascot)
			return ((Mascot) mascotSource).structuredTreeIntervalsInput.get();
		return ((MascotLogPflag) mascotSource).structuredTreeIntervalsInput.get();
	}

	private Dynamics getDynamics() {
		if (mascotSource instanceof Mascot)
			return ((Mascot) mascotSource).dynamicsInput.get();
		return ((MascotLogPflag) mascotSource).dynamicsInput.get();
	}

    @Override
    public void initAndValidate() {
    	if (mascotInput.get() == null && mascotLogPflagInput.get() == null)
    		throw new IllegalArgumentException("Exactly one of mascot and mascotLogPflag must be specified.");
    	if (mascotInput.get() != null && mascotLogPflagInput.get() != null)
    		throw new IllegalArgumentException("Only one of mascot and mascotLogPflag must be specified.");
    	mascotSource = mascotInput.get() != null ? mascotInput.get() : mascotLogPflagInput.get();
    	tree = getStructuredTreeIntervals().treeInput.get();

        if (parameterInput.get().size() == 0 && clockModelInput.get() == null) {
        	someMetaDataNeedsLogging = false;
        	return;
        }
    	someMetaDataNeedsLogging = true;
        if (clockModelInput.get() != null) {
        	substitutions = substitutionsInput.get();
        }

        if (maxStateInput.get() != null){
        	takeMax = maxStateInput.get();

        }

        int dp = decimalPlacesInput.get();

        if (dp < 0) {
            df = null;
        } else {
            df = new DecimalFormat("#."+new String(new char[dp]).replace('\0', '#'));
            df.setRoundingMode(RoundingMode.HALF_UP);
        }

        states = 0;
    }

    @Override
    public void init(PrintStream out) {
    	getStructuredTreeIntervals().treeInput.get().init(out);
    	states = getDynamics().getDimension();
   }

    public void log(int nSample, PrintStream out) {
    	log((long) nSample, out);
    }


    @Override
    public void log(long nSample, PrintStream out) {
    	states = getDynamics().getDimension();

    	try {
			CalculateNodeStates();
		} catch (Exception e) {
			e.printStackTrace();
		}

    	used = new boolean[stateProbabilities.length];
    	report = false;
        List<Function> metadata = parameterInput.get();
        for (int i = 0; i < metadata.size(); i++) {
        	if (metadata.get(i) instanceof StateNode) {
        		metadata.set(i, ((StateNode) metadata.get(i)).getCurrent());
        	}
        }
        BranchRateModel.Base branchRateModel = clockModelInput.get();
        out.print("tree STATE_" + nSample + " = ");
        getStructuredTreeIntervals().treeInput.get().getRoot().sort();
        root = getStructuredTreeIntervals().treeInput.get().getRoot();
        out.print(toNewick(root, metadata, branchRateModel));
        out.print(";");

        for (int i = 0; i < used.length; i++)
        	if(!used[i])
        		System.err.println("not all nodes used");
        if (report)
        	System.err.println("error in node numbers");
    }
    Node root;

	private void appendDouble(StringBuffer buf, double d) {
        if (df == null) {
            buf.append(d);
        } else {
            buf.append(df.format(d));
        }
    }

    String toNewick(Node node, List<Function> metadataList, BranchRateModel.Base branchRateModel) {
        if (maxStateInput.get() != null){
        	takeMax = maxStateInput.get();

        }
        StringBuffer buf = new StringBuffer();
        if (node.getLeft() != null) {
            buf.append("(");
            buf.append(toNewick(node.getLeft(), metadataList, branchRateModel));
            if (node.getRight() != null) {
                buf.append(',');
                buf.append(toNewick(node.getRight(), metadataList, branchRateModel));
            }
            buf.append(")");
        } else {
            buf.append(node.getNr() + 1);
        }
        if (!node.isLeaf()) {
        	if (leftID[node.getNr()-nrSamples] != node.getRight().getNr() && leftID[node.getNr()-nrSamples] != node.getLeft().getNr()){
        		report = true;
        		System.out.println("wrong nr of internal node: "
        				+ leftID[node.getNr()-nrSamples] + " " + node.getLeft().getNr() + " "
        				+ rightID[node.getNr()-nrSamples] + " " + node.getRight().getNr());
        		System.out.println(node.isRoot() + " " + node.getParent().isRoot() + " " + node.getLeft().getID() + " " + node.getRight().getID());
        		System.out.println(getStructuredTreeIntervals().treeInput.get());
        		System.out.println(node.getTree());
        		System.exit(0);
        	}
        	if (!takeMax){
		        buf.append("[&");

		        double[] stateProbs = getStateProb(node.getNr());

		        for (int i = 0 ; i < states; i++)
		        	buf.append(String.format(Locale.US, "%s=%.3f,", getDynamics().getStringStateValue(i), stateProbs[i]));

		        buf.append("max=");
		        buf.append(String.format("%s",
		        		getDynamics().getStringStateValue(whichMax(stateProbs))));

		        if (branchRateModel != null) {
		            buf.append(",rate=");
	                appendDouble(buf, branchRateModel.getRateForBranch(node));
		        }
		        buf.append(']');
        	}else{
		        buf.append("[&max" + type + "=");
		        double[] stateProbs = getStateProb(node.getNr());

		        buf.append(String.format("%d", whichMax(stateProbs) ));

		        if (branchRateModel != null) {
		            buf.append(",rate=");
	                appendDouble(buf, branchRateModel.getRateForBranch(node));
		        }

		        buf.append(']');
        	}
        }else{
			String sampleID = node.getID();
			String[] splits = sampleID.split("_");
			int sampleState;

			if(getDynamics().typeTraitInput.get()!=null){
				sampleState = getDynamics().getValue(node.getID());
			}

			else{
				sampleState = Integer.parseInt(splits[splits.length-1]);
			}
			if (!takeMax){

		        buf.append("[&");

		        for (int i = 0 ; i < states; i++){
		        	if (sampleState != i) buf.append(String.format("%s=0,", getDynamics().getStringStateValue(i)));
		        	if (sampleState == i) buf.append(String.format("%s=1,", getDynamics().getStringStateValue(i)));
		        }
		        buf.append("max=");

		        buf.append(String.format("%s",
		        		getDynamics().getStringStateValue(sampleState)) );
		        buf.append(']');
        	}else{
		        buf.append("[&max" + type + "=");

		        buf.append(String.format("%d", sampleState ));
		        buf.append(']');
        	}
        }

        buf.append(":");
        if (substitutions) {
            appendDouble(buf, node.getLength() * branchRateModel.getRateForBranch(node));
        } else {
            appendDouble(buf, node.getLength());
        }
        return buf.toString();
    }

	@Override
    public void close(PrintStream out) {
		getStructuredTreeIntervals().treeInput.get().close(out);
    }

	//===================================================
	// Calculate the state of all nodes using the up-down
	// algorithm
	//===================================================
	public int samples;
	public int nrSamples;
	public double[][] stateProbabilities;
	public double[][] stateProbabilitiesDown;
	public double[][][] TransitionProbabilities;
	public int[] leftID;
	public int[] rightID;

    public int nrLineages;

    private double[] migrationRates;
    private int[] indicators;
    protected double[] coalescentRates;

    protected ArrayList<Integer> activeLineages;
	private double[] linProbs;
	private double[] transitionProbs;

    private double maxTolerance = 1e-5;
    private boolean recalculateLogP;


    public void CalculateNodeStates() throws Exception{
    	getStructuredTreeIntervals().calculateIntervals();
    	getStructuredTreeIntervals().swap();

    	leftID = new int[getStructuredTreeIntervals().getSampleCount()];
    	rightID = new int[getStructuredTreeIntervals().getSampleCount()];

    	stateProbabilities = new double[getStructuredTreeIntervals().getSampleCount()][];
    	stateProbabilitiesDown = new double[getStructuredTreeIntervals().getSampleCount()][];
    	TransitionProbabilities = new double[getStructuredTreeIntervals().getSampleCount()*2][][];
        nrSamples = getStructuredTreeIntervals().getSampleCount() + 1;

        activeLineages = new ArrayList<Integer>();

        nrLineages = 0;
        linProbs = new double[0];
        transitionProbs = new double[0];

        int treeInterval = 0, ratesInterval = 0;
        double nextEventTime = 0.0;
		coalescentRates = getDynamics().getCoalescentRate(ratesInterval);
        migrationRates = getDynamics().getBackwardsMigration(ratesInterval);
		indicators = getDynamics().getIndicators(ratesInterval);
        double nextTreeEvent = getStructuredTreeIntervals().getInterval(treeInterval);
        double nextRateShift = getDynamics().getInterval(ratesInterval);

        int currTreeInterval = 0;
        do {
        	nextEventTime = Math.min(nextTreeEvent, nextRateShift);

        	if (nextEventTime > 0) {
                if(recalculateLogP){
    				System.err.println("ode calculation stuck, reducing tolerance, new tolerance= " + maxTolerance);
    				maxTolerance *=0.1;
                	CalculateNodeStates();
                	return;
                }
                if(stepSizeInput.get()!=null){
    	        	double[] probs_for_ode = new double[linProbs.length + transitionProbs.length];
    	        	double[] oldLinProbs = new double[linProbs.length + transitionProbs.length];

                    for (int i = 0; i<linProbs.length; i++)
                    	oldLinProbs[i] = linProbs[i];
                    for (int i = linProbs.length; i < transitionProbs.length; i++)
                    	oldLinProbs[i] = transitionProbs[i-linProbs.length];

	                FirstOrderIntegrator integrator = new ClassicalRungeKuttaIntegrator(stepSizeInput.get());
	                FirstOrderDifferentialEquations ode = new MascotODEUpDown(migrationRates, coalescentRates, nrLineages , coalescentRates.length);

	                try {
	                	integrator.integrate(ode, 0, oldLinProbs, nextEventTime, probs_for_ode);
	                }catch(Exception e){
	                	System.out.println(e);
	                	System.out.println("expection");
	                	System.exit(0);
	                	recalculateLogP = true;
	                }

	                for (int i = 0; i<linProbs.length; i++)
	            		linProbs[i] = probs_for_ode[i];
	                for (int i = linProbs.length; i < transitionProbs.length; i++)
	            		transitionProbs[i-linProbs.length] = probs_for_ode[i];
	        	}else {
		        	double[] linProbs_tmp = new double[linProbs.length + transitionProbs.length];
		        	double[] linProbs_tmpdt = new double[linProbs.length + transitionProbs.length];
		        	double[] linProbs_tmpddt = new double[linProbs.length + transitionProbs.length];
		        	double[] linProbs_tmpdddt = new double[linProbs.length + transitionProbs.length];

                    for (int i = 0; i<linProbs.length; i++)
                    	linProbs_tmp[i] = linProbs[i];

                    for (int i = linProbs.length; i < (transitionProbs.length+linProbs.length); i++)
                    	linProbs_tmp[i] = transitionProbs[i-linProbs.length];

                    Euler2ndOrderTransitions euler;
	        		if (getDynamics().hasIndicators)
	        			euler = new Euler2ndOrderTransitions(migrationRates, indicators, coalescentRates, nrLineages , coalescentRates.length, epsilonInput.get(), maxStepInput.get());
	        		else
	        			euler = new Euler2ndOrderTransitions(migrationRates, coalescentRates, nrLineages , coalescentRates.length, epsilonInput.get(), maxStepInput.get());

		        	linProbs[linProbs.length-1] = 0;
		        	euler.calculateValues(nextEventTime, linProbs_tmp, linProbs_tmpdt, linProbs_tmpddt, linProbs_tmpdddt);

	                for (int i = 0; i<linProbs.length; i++)
	            		linProbs[i] = linProbs_tmp[i];
	                for (int i = linProbs.length; i < linProbs_tmp.length; i++)
	            		transitionProbs[i-linProbs.length] = linProbs_tmp[i];

	        	}
			}

        	if (nextTreeEvent <= nextRateShift){
 	        	if (getStructuredTreeIntervals().getIntervalType(treeInterval) == IntervalType.COALESCENT) {
 	        		nrLineages--;
	        		normalizeLineages();
        			coalesce(treeInterval);
	        	}

 	       		if (getStructuredTreeIntervals().getIntervalType(treeInterval) == IntervalType.SAMPLE) {
 	       			nrLineages++;
 	       			if (linProbs.length > 0)
 	       				normalizeLineages();
 	       			sample(treeInterval);
	       		}

 	       		treeInterval++;
        		nextRateShift -= nextTreeEvent;
        		try{
        			nextTreeEvent = getStructuredTreeIntervals().getInterval(treeInterval);
        		}catch(Exception e){
        			break;
        		}
        	}else{
        		ratesInterval++;
        		coalescentRates = getDynamics().getCoalescentRate(ratesInterval);
                migrationRates = getDynamics().getBackwardsMigration(ratesInterval);
        		indicators = getDynamics().getIndicators(ratesInterval);
        		nextTreeEvent -= nextRateShift;
 	       		nextRateShift = getDynamics().getInterval(ratesInterval);
        	}

        }while(nextTreeEvent <= Double.POSITIVE_INFINITY);
        currTreeInterval = getStructuredTreeIntervals().getIntervalCount()-1;

        do{
		  	if (getStructuredTreeIntervals().getIntervalType(currTreeInterval) == IntervalType.COALESCENT) {
		  		coalesceDown(currTreeInterval);
		   	}
		  	currTreeInterval--;
        }while(currTreeInterval>=0);
    }


    private double normalizeLineages(){
    	if (linProbs==null)
    		return 0.0;

    	double interval = 0.0;
    	for (int i = 0; i < linProbs.length/states; i++){
    		double lineProbs = 0.0;
    		for (int j = 0; j < states; j++)
    			if (linProbs[i*states+j]>=0.0){
    				lineProbs += linProbs[i*states+j];
    			}else{
    				recalculateLogP = true;
    				return Math.log(1.0);
    			}
    		for (int j = 0; j < states; j++){
    			linProbs[i*states+j] = linProbs[i*states+j]/lineProbs;
    		}
    		interval +=lineProbs;
    	}

		return Math.log(interval/(linProbs.length/states));

    }

    private void sample(int currTreeInterval) {
		int incomingLines = getStructuredTreeIntervals().getLineagesAdded(currTreeInterval);
		int newLengthLineages = linProbs.length + 1*states;
		int newLengthTransitions = transitionProbs.length + 1*states*states;

		double[] linProbsNew = new double[newLengthLineages];
		double[] transitionProbsNew = new double[newLengthTransitions];

		for (int i = 0; i < linProbs.length; i++)
			linProbsNew[i] = linProbs[i];

		for (int i = 0; i < transitionProbs.length; i++)
			transitionProbsNew[i] = transitionProbs[i];

		int currPositionLineages = linProbs.length;
		int currPositionTransitions = transitionProbs.length;

		if (getDynamics().typeTraitInput.get()!=null){
			int  l = incomingLines; {
				activeLineages.add(l);
				int sampleState = (int) getDynamics().getValue(tree.getNode(l).getID());
				for (int i = 0; i< states; i++){
					if (i == sampleState){
						linProbsNew[currPositionLineages] = 1.0;currPositionLineages++;
					}
					else{
						linProbsNew[currPositionLineages] = 0.0;currPositionLineages++;
					}
				}
				for (int s = 0; s < states; s++){
					for (int i = 0; i < states; i++){
						if (i == s){
							transitionProbsNew[currPositionTransitions] = 1.0;
							currPositionTransitions++;
						}else{
							transitionProbsNew[currPositionTransitions] = 0.0;
							currPositionTransitions++;
						}
					}
				}
			}
		}else{
			int l = incomingLines; {
				activeLineages.add(l);
				String sampleID = tree.getNode(l).getID();
				int sampleState = 0;
				if (states > 1){
					String[] splits = sampleID.split("_");
					sampleState = Integer.parseInt(splits[splits.length-1]);
				}
				for (int i = 0; i< states; i++){
					if (i == sampleState){
						linProbsNew[currPositionLineages] = 1.0;currPositionLineages++;
					}
					else{
						linProbsNew[currPositionLineages] = 0.0;currPositionLineages++;
					}
				}
				for (int s = 0; s < states; s++){
					for (int i = 0; i < states; i++){
						if (i == s){
							transitionProbsNew[currPositionTransitions] = 1.0;
							currPositionTransitions++;
						}else{
							transitionProbsNew[currPositionTransitions] = 0.0;
							currPositionTransitions++;
						}
					}
				}
			}
		}
		linProbs = linProbsNew;
		transitionProbs = transitionProbsNew;
    }

    private void coalesce(int currTreeInterval) {
    	for (int i = 0; i < nrLineages*states; i++){
    		double lineProbs = 0.0;
    		for (int j = 0; j < states; j++){
    			if (transitionProbs.length>=0.0){
    				lineProbs += transitionProbs[i*states+j];
    			}else{
    				System.err.println("transition probability smaller than 0 (or NaN before normalizing)");
    				System.exit(0);
    			}
    		}
    		for (int j = 0; j < states; j++)
    			transitionProbs[i*states+j] = transitionProbs[states*i+j]/lineProbs;
    	}

    	int [] coalLines = new int[] {
    			getStructuredTreeIntervals().getLineagesRemoved(currTreeInterval,0),
    			getStructuredTreeIntervals().getLineagesRemoved(currTreeInterval,1)
    	};
    	final int daughterIndex1 = activeLineages.indexOf(coalLines[0]);
		final int daughterIndex2 = activeLineages.indexOf(coalLines[1]);
		if (daughterIndex1 == -1 || daughterIndex2 == -1) {
			System.out.println("daughter lineages at coalescent event not found");
    		System.exit(0);
		}
		double[] lambda = new double[states];
		double lambdaSum = 0;

        for (int k = 0; k < states; k++) {
        	Double pairCoalRate = coalescentRates[k] * linProbs[daughterIndex1*states + k] * linProbs[daughterIndex2*states + k];
			if (!Double.isNaN(pairCoalRate)){
				lambda[k] =  pairCoalRate;
				lambdaSum += pairCoalRate;
			}
        }

        activeLineages.add(tree.getNode(coalLines[0]).getParent().getNr());

		double[] pVec = new double[states];
		for (int i = 0; i < pVec.length; i++)
			pVec[i] = lambda[i]/lambdaSum;

		stateProbabilities[tree.getNode(coalLines[0]).getParent().getNr() - nrSamples] = pVec;

		double[][] tP1 = new double[states][states];
		for (int i = 0; i< states; i++){
			for (int j = 0; j< states; j++){
				tP1[i][j] = transitionProbs[daughterIndex1*states*states+i*states+j];
			}
		}

		double[][] tP2 = new double[states][states];
		for (int i = 0; i< states; i++){
			for (int j = 0; j< states; j++){
				tP2[i][j] = transitionProbs[daughterIndex2*states*states+i*states+j];
			}
		}

		double[] linProbsNew  = new double[linProbs.length - states];

		int linCount = 0;
		for (int i = 0; i <= nrLineages; i++){
			if (i != daughterIndex1 && i != daughterIndex2){
				for (int j = 0; j < states; j++){
					linProbsNew[linCount*states + j] = linProbs[i*states + j];
				}
				linCount++;
			}
		}
		for (int j = 0; j < states; j++){
			linProbsNew[linCount*states + j] = pVec[j];
		}
		linProbs = linProbsNew;

		double[] transitionProbsNew  = new double[transitionProbs.length - states*states];

		linCount = 0;
		for (int i = 0; i <= nrLineages; i++){
			if (i != daughterIndex1 && i != daughterIndex2){
				for (int j = 0; j < states; j++)
					for (int k = 0; k < states; k++)
						transitionProbsNew[linCount*states*states+j*states+k]
								= transitionProbs[i*states*states+j*states+k];
				linCount++;
			}
		}

		for (int j = 0; j < states; j++)
			for (int k = 0; k < states; k++)
				if (j==k)
					transitionProbsNew[linCount*states*states+j*states+k] = 1.0;
				else
					transitionProbsNew[linCount*states*states+j*states+k] = 0.0;

		transitionProbs = transitionProbsNew;

		TransitionProbabilities[coalLines[0]] = tP1;
		TransitionProbabilities[coalLines[1]] = tP2;

		if (daughterIndex1>daughterIndex2){
			activeLineages.remove(daughterIndex1);
			activeLineages.remove(daughterIndex2);
		}else{
			activeLineages.remove(daughterIndex2);
			activeLineages.remove(daughterIndex1);
		}

		if(tree.getNode(coalLines[0]).getParent().getNr() != getStructuredTreeIntervals().treeInput.get().getNode(coalLines[1]).getParent().getNr())
			System.err.println("wrong daughter parent");
		if(tree.getNode(coalLines[1]).getParent().getNr() != getStructuredTreeIntervals().treeInput.get().getNode(coalLines[0]).getParent().getNr())
			System.err.println("wrong daughter parent");
		if(tree.getNode(coalLines[1]).getParent().getNr() != tree.getNode(coalLines[0]).getParent().getNr())
			System.err.println("coalescent nodes don't have the same parent");

		leftID[tree.getNode(coalLines[0]).getParent().getNr() - nrSamples] = tree.getNode(coalLines[0]).getParent().getLeft().getNr();
		rightID[tree.getNode(coalLines[0]).getParent().getNr() - nrSamples] = tree.getNode(coalLines[0]).getParent().getRight().getNr();
    }

    private void coalesceDown(int currTreeInterval) {
		int parentLines = getStructuredTreeIntervals().getLineagesAdded(currTreeInterval);
		Node parentNode = tree.getNode(parentLines);

		if (!parentNode.isRoot()){
			double[] start = stateProbabilities[parentNode.getNr() - nrSamples];
			double[] end = stateProbabilitiesDown[parentNode.getParent().getNr() - nrSamples];
			double[][] flow = TransitionProbabilities[parentNode.getNr()];
			double[] otherSideInfo = new double[states];
			for (int a = 0; a < states; a++) {
				double sum = 0;
				for (int b = 0; b < states; b++) {
					sum += start[b] * flow[b][a];
				}
				otherSideInfo[a] = end[a]/sum;
				if (Double.isNaN(otherSideInfo[a]))
					otherSideInfo[a] = 0;
			}

			double[] conditional = new double[states];
			double condsum = 0;
			for (int a = 0; a < states; a++) {
				double sum = 0;
				for (int b = 0; b < states; b++) {
					sum += flow[a][b] * otherSideInfo[b];
				}
				conditional[a] = sum *start[a];
				condsum += conditional[a];
			}
			for (int a = 0; a < states; a++)
				conditional[a] /= condsum;

			stateProbabilitiesDown[parentNode.getNr() - nrSamples] = conditional;
		}else{
			stateProbabilitiesDown[parentNode.getNr() - nrSamples] = stateProbabilities[parentNode.getNr() - nrSamples];
    	}
	}


	public double[] getStateProb(int nr) {
		if(useUpDown.get()){
			used[nr - nrSamples] = true;
			return stateProbabilitiesDown[nr - nrSamples] ;
		}else{
			used[nr - nrSamples] = true;
			return stateProbabilities[nr - nrSamples] ;
		}
	}

	public double[] getStateProbOnly(int nr) {
		if(useUpDown.get()){
			return stateProbabilitiesDown[nr - nrSamples] ;
		}else{
			return stateProbabilities[nr - nrSamples] ;
		}
	}

	public int whichMax(double[] stateProbs) {
		double max_val = -1;
		int max_ind = 1;
		for (int i = 0; i < stateProbs.length;i++) {
			if (stateProbs[i]>max_val) {
				max_val = stateProbs[i];
				max_ind = i;
			}
		}
		return max_ind;
	}

    public void calcForTest() {
    	states = getDynamics().getDimension();
    	try {
			CalculateNodeStates();
		} catch (Exception e) {
			e.printStackTrace();
		}
    }

}
