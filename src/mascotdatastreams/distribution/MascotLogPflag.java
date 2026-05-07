package mascotdatastreams.distribution;


import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.IntervalType;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.CalculationNode;
import mascot.distribution.MascotNative2;
import mascot.distribution.StructuredTreeDistribution;
import mascot.distribution.StructuredTreeIntervals;
import mascot.dynamics.Dynamics;
import mascot.ode.*;

import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * @author Nicola Felix Mueller
 */

@Description("Calculates the probability of a beast.tree using under the framework of Mueller (2017).")
@Citation("Nicola F. Müller, David A. Rasmussen, Tanja Stadler (2017)\n  The Structured Coalescent and its Approximations.\n  Mol Biol Evol 2017 msx186. doi: 10.1093/molbev/msx186")
public class MascotLogPflag extends StructuredTreeDistribution {

	public static boolean debug = false;

	/**
	 * Mechanism-level metric for benchmarking: how many ODE steps were taken
	 * across all calculateLogP / calculateLogP_maxInterval calls since the
	 * counter was last reset. Incremented in both doEuler and doEulerAtTime.
	 * Process-wide static — only meaningful when running one likelihood at a
	 * time (true for our benchmark harnesses).
	 */
	public static long doEulerCallCount = 0;
	public static long getAndResetDoEulerCallCount() {
		long v = doEulerCallCount;
		doEulerCallCount = 0;
		return v;
	}

	public Input<Dynamics> dynamicsInput = new Input<>("dynamics", "Input of rates", Input.Validate.REQUIRED);
	public Input<Double> epsilonInput = new Input<>("epsilon", "step size for the RK4 integration",0.001);
	public Input<Double> maxStepInput = new Input<>("maxStep", "step size for the RK4 integration", Double.POSITIVE_INFINITY);

	public Input<Boolean> cacheInput = new Input<>("useCache", "use cache to speed things up (may be fragile)", false);
	public Input<Boolean> computeLikelihoodInput = new Input<>("compute_likelihood", "If true, compute the Mascot tree likelihood; if false, return 0 after initialization.", true);

	public Input<Boolean> useMaxIntervalInput = new Input<>("useMaxInterval",
			"If true, ignore the dynamics grid for ODE subdivision; subdivide each tree-event interval into ceil(L / maxInterval) equal pieces and evaluate rates at each subinterval midpoint. Requires implementation=java, useCache=false, and StructuredSkylinePrevalence dynamics.",
			false);
	public Input<Double> maxIntervalInput = new Input<>("maxInterval",
			"Maximum length of any single ODE step when useMaxInterval=true.",
			Double.POSITIVE_INFINITY);

	enum MascotImplementation {java, indicators, allnative};
	public Input<MascotImplementation> implementationInput = new Input<>("implementation", "implementation, one of " + MascotImplementation.values().toString(),
			MascotImplementation.allnative, MascotImplementation.values());

	public int samples;
	public int nrSamples;
	public double[][] stateProbabilities;

	protected int nrLineages;

    // current rates
    //private double[] migrationRates;
    //private int [] indicators;
    protected double[] coalescentRates;


    // Set up for lineage state probabilities
    protected ArrayList<Integer> activeLineages;
    protected double[] linProbs;
    protected double[] linProbsNew;
    protected int linProbsLength;
    protected int states;

    // store the linProbs, multiplicators and logP's at coalescent points in jagged arrays from last time
    protected double[] coalLinProbs;
    protected int [] coalLinProbsLengths;
    protected double[] coalLogP;
    protected int[] coalRatesInterval;
    //private ArrayList<ArrayList<Integer>> coalActiveLineages;

    // deep store the things above for MCMC
    protected double[] storeLinProbs;
    protected int [] storedCoalLinProbsLengths;
    protected double[] storeLogP;
    protected int[] storeRatesInterval;
    //private ArrayList<ArrayList<Integer>> storeActiveLineages;


    protected double [] nextTreeEvents;
    protected double [] storedNextTreeEvents;
    protected double [] nextRateShifts;
    protected double [] storedNextRateShifts;
	//private int [] treeIntervalNrs;
	//private int [] storedTreeIntervalNrs;
	//private int [] lineagesAddded;
	//private int [] storedLineagesAddded;

    // check if this is the first calculation
    protected int first = 0;


	// maximum integration error tolerance
    protected double maxTolerance = 1e-3;
    protected boolean recalculateLogP;
	Euler2ndOrderBase euler;
	TreeInterface tree;
	public Dynamics dynamics;
	StructuredTreeIntervals treeIntervals;

	int [] nodeType;

	MascotNative2 mascotImpl = null;
	boolean useCache;
	boolean computeLikelihood;

    @Override
    public void initAndValidate(){
    	useCache = cacheInput.get();
    	computeLikelihood = computeLikelihoodInput.get();
    	dynamics = dynamicsInput.get();
    	treeIntervals = structuredTreeIntervalsInput.get();
//    	tree = treeInput.get();
    	if (tree == null) {
    		tree = treeIntervals.treeInput.get();
    	}
    	treeIntervals.calculateIntervals();
    	stateProbabilities = new double[treeIntervals.getSampleCount()][states];
        nrSamples = treeIntervals.getSampleCount() + 1;
        states = dynamics.getDimension();

    	int intCount = treeIntervals.getIntervalCount();

    	// initialize storing arrays and ArrayLists
    	coalLinProbs = new double[intCount * intCount * states];
    	storeLinProbs = new double[intCount * intCount * states];
    	coalLinProbsLengths = new int[intCount];
    	storedCoalLinProbsLengths = new int[intCount];
    	coalLogP = new double[intCount];
    	storeLogP = new double[intCount];
    	coalRatesInterval = new int[intCount];
    	storeRatesInterval = new int[intCount];
    	//coalActiveLineages = new ArrayList<>();
    	nextTreeEvents = new double[intCount];
    	nextRateShifts = new double[intCount];
    	storedNextTreeEvents = new double[intCount];
    	storedNextRateShifts = new double[intCount];
    	parents = new int[intCount];
    	//treeIntervalNrs = new int[intCount];
    	//storedTreeIntervalNrs = new int[intCount];
    	//lineagesAddded = new int[intCount];
    	//storedLineagesAddded = new int[intCount];

    	//ArrayList<Integer> emptyList = new ArrayList<>();
    	//for (int i = 0; i <= intCount; i++) coalActiveLineages.add(emptyList);

    	activeLineages = new ArrayList<>();

    	int MAX_SIZE = intCount * states;
    	linProbs_for_ode = new double[MAX_SIZE];
    	linProbs_tmp = new double[MAX_SIZE];
    	linProbs = new double[MAX_SIZE];
    	linProbsNew = new double[MAX_SIZE];

		nodeType = new int[tree.getNodeCount()];
    	if (dynamics.typeTraitInput.get() != null) {
    		nodeType = new int[tree.getNodeCount()];
    		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
    			nodeType[i] = dynamics.getValue(tree.getNode(i).getID());
				if (nodeType[i] >= dynamics.getDimension()) {
					System.err.println("sample discovered for node id='"+tree.getNode(i).getID()+"' with higher state than dimension");
//					System.exit(1);
				}
    		}
    	} else {
    		// TODO: fill in nodeType another way
    	}

    	MascotImplementation imp = implementationInput.get();
    	switch (imp) {
    	case allnative: if (Euler2ndOrderNative.loadLibrary()) {
    		mascotImpl = new MascotNative2(treeIntervals, nodeType, states,epsilonInput.get(), maxStepInput.get(), useCache);
    		break;
    	}
    	case indicators: if (Euler2ndOrderNative.loadLibrary()) {
    		euler = new Euler2ndOrderNative();
        	euler.setup(MAX_SIZE, states, epsilonInput.get(), maxStepInput.get());
        	Log.warning("Using " + euler.getClass().getSimpleName());
    		break;
    	}
    	case java:
    		switch (states) {
    		case 2: euler = new Euler2ndOrder2(); break;
    		case 3: euler = new Euler2ndOrder3(); break;
    		case 4: euler = new Euler2ndOrder4(); break;
    		case 5: euler = new Euler2ndOrder5(); break;
    		case 6: euler = new Euler2ndOrder6(); break;
    		case 7: euler = new Euler2ndOrder7(); break;
    		case 8: euler = new Euler2ndOrder8(); break;
    		case 9: euler = new Euler2ndOrder9(); break;
    		case 10: euler = new Euler2ndOrder10(); break;
    		case 11: euler = new Euler2ndOrder11(); break;
    		case 12: euler = new Euler2ndOrder12(); break;
    		case 13: euler = new Euler2ndOrder13(); break;
    		case 14: euler = new Euler2ndOrder14(); break;
    		case 15: euler = new Euler2ndOrder15(); break;
    		case 16: euler = new Euler2ndOrder16(); break;
    		case 17: euler = new Euler2ndOrder17(); break;
    		case 18: euler = new Euler2ndOrder18(); break;
    		case 19: euler = new Euler2ndOrder19(); break;
    		case 20: euler = new Euler2ndOrder20(); break;
    		case 21: euler = new Euler2ndOrder21(); break;
    		case 22: euler = new Euler2ndOrder22(); break;
    		case 23: euler = new Euler2ndOrder23(); break;
    		case 24: euler = new Euler2ndOrder24(); break;
    		case 25: euler = new Euler2ndOrder25(); break;
    		case 26: euler = new Euler2ndOrder26(); break;
    		case 27: euler = new Euler2ndOrder27(); break;
    		case 28: euler = new Euler2ndOrder28(); break;
    		case 29: euler = new Euler2ndOrder29(); break;
    		case 30: euler = new Euler2ndOrder30(); break;
    		default: euler = new Euler2ndOrder(); break;
    		}

        	euler.setup(MAX_SIZE, states, epsilonInput.get(), maxStepInput.get());
        	Log.warning("Using " + euler.getClass().getSimpleName());
    	}

    	if (useMaxIntervalInput.get()) {
    		if (cacheInput.get())
    			throw new IllegalArgumentException("useMaxInterval=true is incompatible with useCache=true");
    		if (implementationInput.get() != MascotImplementation.java)
    			throw new IllegalArgumentException("useMaxInterval=true requires implementation=\"java\"");
    		if (!(dynamics instanceof StructuredSkylinePrevalence))
    			throw new IllegalArgumentException("useMaxInterval=true requires dynamics of type StructuredSkylinePrevalence");
    		if (maxIntervalInput.get() <= 0.0)
    			throw new IllegalArgumentException("maxInterval must be > 0");

    		// Swap in the vendored Euler classes (mascotdatastreams.ode) which
    		// fix initWithIndicators (upstream forgets to allocate sumDotStates).
    		// doEulerAtTime relies on initWithIndicators to enable the
    		// indicator-based sparsity optimisation in the inner loop.
    		switch (states) {
    		case 2: euler = new mascotdatastreams.ode.Euler2ndOrder2(); break;
    		case 3: euler = new mascotdatastreams.ode.Euler2ndOrder3(); break;
    		case 4: euler = new mascotdatastreams.ode.Euler2ndOrder4(); break;
    		default: euler = new mascotdatastreams.ode.Euler2ndOrder(); break;
    		}
    		euler.setup(MAX_SIZE, states, epsilonInput.get(), maxStepInput.get());
    		Log.warning("MascotLogPflag: useMaxInterval=true, maxInterval=" + maxIntervalInput.get()
    				+ ", euler=" + euler.getClass().getName());
    	}

    }

    double [] linProbs_for_ode;
    double [] linProbs_tmp;
    int [] parents;

    /** Calls package-private MascotNative2.setUpDynamics(Dynamics) via reflection. */
    private static void callMascotNative2SetUpDynamics(MascotNative2 impl, Dynamics dynamics) {
        try {
            Method m = MascotNative2.class.getDeclaredMethod("setUpDynamics", Dynamics.class);
            m.setAccessible(true);
            m.invoke(impl, dynamics);
        } catch (ReflectiveOperationException e) {
            throw new RuntimeException("Failed to call MascotNative2.setUpDynamics (package-private)", e);
        }
    }

    /** Reads StructuredTreeIntervals.firstDirtyInterval (package-private) via reflection. */
    private static int getStiFirstDirtyInterval(StructuredTreeIntervals ti) {
        return getStiFieldInt(ti, "firstDirtyInterval");
    }

    /** Reads StructuredTreeIntervals.lineagesAdded (protected, not visible across packages) via reflection. */
    private static int[] getStiLineagesAdded(StructuredTreeIntervals ti) {
        return getStiFieldIntArray(ti, "lineagesAdded");
    }

    /** Reads StructuredTreeIntervals.lineagesRemoved via reflection. */
    private static int[] getStiLineagesRemoved(StructuredTreeIntervals ti) {
        return getStiFieldIntArray(ti, "lineagesRemoved");
    }

    /** Reads StructuredTreeIntervals.intervals via reflection. */
    private static double[] getStiIntervals(StructuredTreeIntervals ti) {
        try {
            Field f = StructuredTreeIntervals.class.getDeclaredField("intervals");
            f.setAccessible(true);
            return (double[]) f.get(ti);
        } catch (ReflectiveOperationException e) {
            throw new RuntimeException("Failed to read StructuredTreeIntervals.intervals", e);
        }
    }

    /** Reads StructuredTreeIntervals.storedLineagesAdded via reflection. */
    private static int[] getStiStoredLineagesAdded(StructuredTreeIntervals ti) {
        return getStiFieldIntArray(ti, "storedLineagesAdded");
    }

    private static int getStiFieldInt(StructuredTreeIntervals ti, String fieldName) {
        try {
            Field f = StructuredTreeIntervals.class.getDeclaredField(fieldName);
            f.setAccessible(true);
            return f.getInt(ti);
        } catch (ReflectiveOperationException e) {
            throw new RuntimeException("Failed to read StructuredTreeIntervals." + fieldName, e);
        }
    }

    private static int[] getStiFieldIntArray(StructuredTreeIntervals ti, String fieldName) {
        try {
            Field f = StructuredTreeIntervals.class.getDeclaredField(fieldName);
            f.setAccessible(true);
            return (int[]) f.get(ti);
        } catch (ReflectiveOperationException e) {
            throw new RuntimeException("Failed to read StructuredTreeIntervals." + fieldName, e);
        }
    }

    public double calculateLogP() {
    	// newly calculate tree intervals (already done by swap() below)
    	treeIntervals.calculateIntervals();
    	// correctly calculate the daughter nodes at coalescent intervals in the case of
    	// bifurcation or in case two nodes are at the same height
    	treeIntervals.swap();

    	if (useMaxIntervalInput.get()) {
    		return calculateLogP_maxInterval();
    	}

    	if (mascotImpl != null) {
    		Node [] nodes = tree.getNodesAsArray();
    		for (int i = 0; i < nodes.length - 1; i++) {
    			parents[i] = nodes[i].getParent().getNr();
    		}
            if (first == 0 || !dynamics.areDynamicsKnown()) {
            	callMascotNative2SetUpDynamics(mascotImpl, dynamics);
            }
    		if (!computeLikelihood) {
    			first++;
    			return 0;
    		}
    		logP = mascotImpl.calculateLogP(dynamics.isDirtyCalculation() || first == 0,
    				getStiFirstDirtyInterval(treeIntervals),
    				getStiLineagesAdded(treeIntervals),
    				getStiLineagesRemoved(treeIntervals),
    				getStiIntervals(treeIntervals),
    				parents
    				);
    		//System.out.println(logP);
    		first++;
    		return logP;
    	}
        // Set up ArrayLists for the indices of active lineages and the lineage state probabilities
        activeLineages.clear();
        logP = 0;
        nrLineages = 0;
        //linProbs = new double[0];// initialize the tree and rates interval counter
        linProbsLength = 0;
        int treeInterval = 0, ratesInterval = 0;
     	double nextEventTime = 0.0;

        // Time to the next rate shift or event on the tree
        double nextTreeEvent = treeIntervals.getInterval(treeInterval);
        double nextRateShift = dynamics.getInterval(ratesInterval);

        //System.err.println("first = " + first);
        if (first == 270) {
        	//debug = true;
        }
        if (debug) {
        	Log.info.println("##" + getStiFirstDirtyInterval(treeIntervals));
        	Log.info.println("##" + Arrays.toString(getStiLineagesAdded(treeIntervals)));
        	Log.info.println("##" + Arrays.toString(getStiLineagesRemoved(treeIntervals)));
        	Log.info.println("##" + Arrays.toString(getStiIntervals(treeIntervals)));
        }
        if (first == 0 || !dynamics.areDynamicsKnown()) {
        	setUpDynamics();
        }
        if (useCache && first > 0 && !dynamics.isDirtyCalculation()
					&& getStiFirstDirtyInterval(treeIntervals) > 2) {
			// restore the likelihood to last known good place
			int pos0 = -1, pos1 = -1;
			do {

				nextEventTime = nextTreeEvent;

				// Check if the last interval was reached
				if (treeInterval == treeIntervals.getIntervalCount()){
					// Log.warning("Restoring to the finish!");
					logP = coalLogP[coalLogP.length-1];
					return logP;
				}
				boolean isDirty;
				if (treeIntervals.getCoalescentEvents(treeInterval) > 0) { // == IntervalType.COALESCENT) {
					Integer coalLines0 = treeIntervals.getLineagesRemoved(treeInterval, 0);
					Integer coalLines1 = treeIntervals.getLineagesRemoved(treeInterval, 1);
					pos0 = activeLineages.indexOf(coalLines0);
					pos1 = activeLineages.indexOf(coalLines1);
					if (pos0 < 0 || pos1 < 0) {
						System.out.println(coalLines0/*.getNr()*/ + " " + coalLines1/*.getNr()*/ + " " + activeLineages);
						System.out.println("daughter lineages at coalescent event not found");
	// 	           		return Double.NEGATIVE_INFINITY;
						throw new RuntimeException("coalesceX went wrong at 1");
					}
					if (pos0 > pos1) {
						activeLineages.remove(pos0);
						activeLineages.remove(pos1);
					} else {
						activeLineages.remove(pos1);
						activeLineages.remove(pos0);
					}
					int newLineage = tree.getNode(coalLines0).getParent().getNr();
					activeLineages.add(newLineage);
					isDirty = getStiStoredLineagesAdded(treeIntervals)[treeInterval] != newLineage;
				} else { // (treeIntervals.getIntervalType(treeInterval) == IntervalType.SAMPLE) {
					int incomingLines = treeIntervals.getLineagesAdded(treeInterval);
					activeLineages.add(incomingLines);
					isDirty = getStiStoredLineagesAdded(treeIntervals)[treeInterval] != incomingLines;
				}


				if (isDirty || treeInterval+1 == getStiFirstDirtyInterval(treeIntervals)) {
	//    				treeIntervals.intervalIsDirty(treeInterval+1) || isDirty) {
					if (treeInterval <= 2) {
						treeInterval = 0;
						ratesInterval = 0;
						activeLineages.clear();
						logP = 0;
						nextEventTime = 0.0;
						nextTreeEvent = treeIntervals.getInterval(treeInterval);
						nextRateShift = dynamics.getInterval(ratesInterval);
						break;
					}

	// 	    			Log.warning.print("Restore I ");
					activeLineages.remove(activeLineages.size() - 1);
					if (treeIntervals.getCoalescentEvents(treeInterval) > 0) { // == IntervalType.COALESCENT) {
						Integer coalLines0 = treeIntervals.getLineagesRemoved(treeInterval, 0);
						Integer coalLines1 = treeIntervals.getLineagesRemoved(treeInterval, 1);
						if (pos0 > pos1) {
							activeLineages.add(pos1, coalLines1);
							activeLineages.add(pos0, coalLines0);
						} else {
							activeLineages.add(pos0, coalLines0);
							activeLineages.add(pos1, coalLines1);
						}
					}

					treeInterval++;
					ratesInterval = restoreNode(treeInterval-2);
					nextTreeEvent = nextTreeEvents[treeInterval-1];
					nextRateShift = nextRateShifts[treeInterval-1];
					treeInterval--;
					break;
				}
				treeInterval++;
				nextRateShift -= nextTreeEvent;
				try{
					nextTreeEvent = treeIntervals.getInterval(treeInterval);
				} catch(Exception e) {
					break;
				}
	//        	} else {
	//        		ratesInterval++;
	//        		nextTreeEvent -= nextRateShift;
	// 	       		nextRateShift = dynamics.getInterval(ratesInterval);
	//        	}
			} while(nextTreeEvent <= Double.POSITIVE_INFINITY);
		}


		coalescentRates = dynamics.getCoalescentRate(ratesInterval);
        //migrationRates = dynamics.getBackwardsMigration(ratesInterval);
		//indicators = dynamics.getIndicators(ratesInterval);
		nrLineages = activeLineages.size();
		linProbsLength = nrLineages * states;

		if (!computeLikelihood) {
			first++;
			return 0;
		}

        // Calculate the likelihood
        do {
        	nextEventTime = Math.min(nextTreeEvent, nextRateShift);
        	if (nextEventTime > 0) {													// if true, calculate the interval contribution
                if (recalculateLogP) {
    				System.err.println("ode calculation stuck, reducing tolerance, new tolerance= " + maxTolerance);
    				maxTolerance *=0.9;
    		    	recalculateLogP = false;
    				System.exit(0);
                	return calculateLogP();
                }
        		logP += doEuler(nextEventTime, ratesInterval);
        	}

        	if (nextTreeEvent <= nextRateShift){
 	        	if (treeIntervals.getIntervalType(treeInterval) == IntervalType.COALESCENT) {
// 	        		System.out.print(String.format("%.3f ", nextTreeEvent));
//	        		logP += normalizeLineages(linProbs);									// normalize all lineages before event
 	        		nrLineages--;													// coalescent event reduces the number of lineages by one
	        		logP += coalesce(treeInterval, ratesInterval, nextTreeEvent, nextRateShift);	  				// calculate the likelihood of the coalescent event
	        	}
 	       		
 	       		if (treeIntervals.getIntervalType(treeInterval) == IntervalType.SAMPLE) {
 	       			//if (linProbsLength > 0)
 	       			//	logP += normalizeLineages(linProbs);								// normalize all lineages before event
	       			nrLineages++;													// sampling event increases the number of lineages by one
 	       			sample(treeInterval, ratesInterval, nextTreeEvent, nextRateShift);							// calculate the likelihood of the sampling event if sampling rate is given
	       		}

 	       		treeInterval++;
        		nextRateShift -= nextTreeEvent;
        		try{
        			nextTreeEvent = treeIntervals.getInterval(treeInterval);
        		}catch(Exception e){
        			break;
        		}
        	} else {
        		ratesInterval++;
        		coalescentRates = dynamics.getCoalescentRate(ratesInterval);
                //migrationRates = dynamics.getBackwardsMigration(ratesInterval);
        		//indicators = dynamics.getIndicators(ratesInterval);
        		nextTreeEvent -= nextRateShift;
 	       		nextRateShift = dynamics.getInterval(ratesInterval);
        	}
        	if (logP == Double.NEGATIVE_INFINITY) {
        		return logP;
        	}
            if (debug) {
            	Log.info(treeInterval + " " + ratesInterval + " " + logP);
            }
        } while(nextTreeEvent <= Double.POSITIVE_INFINITY);

        first++;
//        System.err.println(logP);
		return logP;
    }

    // Per-mascotshifts-cell rate cache for calculateLogP_maxInterval. Mirrors
    // what setUpDynamics builds for the legacy path: rates are evaluated once
    // per epoch midpoint (per likelihood call when dynamics are dirty), then
    // looked up by cell index inside the inner loop. Avoids the per-step
    // double[] allocation and the spline binary-search that the original
    // getCoalescentRateAtTime / getBackwardsMigrationAtTime do.
    private double[][] mxCoalRates;
    private double[][] mxMigRates;
    private double[]   mxRateShiftEnds;   // cumulative time at the END of each cell
    private int[]      mxIndicators;      // time-invariant indicators
    private int        mxNumCells;

    /** Recompute the per-cell rate cache. Called when dynamics are dirty. */
    private void precomputeMaxIntervalCache() {
        int n = dynamics.getEpochCount();
        if (mxCoalRates == null || mxCoalRates.length != n) {
            mxCoalRates = new double[n][];
            mxMigRates = new double[n][];
            mxRateShiftEnds = new double[n];
        }
        double cum = 0.0;
        for (int i = 0; i < n; i++) {
            mxCoalRates[i] = dynamics.getCoalescentRate(i);
            mxMigRates[i] = dynamics.getBackwardsMigration(i);
            cum += dynamics.getInterval(i);
            mxRateShiftEnds[i] = cum;
        }
        mxIndicators = dynamics.getIndicators(0);
        mxNumCells = n;
    }

    /**
     * Max-interval variant of calculateLogP. Subdivides each tree-event interval
     * into ceil(L / maxInterval) equal sub-intervals; rates inside each
     * sub-interval are taken from the precomputed per-cell cache (constant
     * across each mascotshifts cell, matching the legacy path's resolution).
     * Caching is disabled; only the java implementation path and
     * StructuredSkylinePrevalence dynamics are supported (validated in
     * initAndValidate).
     */
    public double calculateLogP_maxInterval() {
        double maxInt = maxIntervalInput.get();

        activeLineages.clear();
        logP = 0;
        nrLineages = 0;
        linProbsLength = 0;

        int treeInterval = 0;
        double currentTime = 0.0;
        double nextTreeEvent;
        try {
            nextTreeEvent = treeIntervals.getInterval(treeInterval);
        } catch (Exception e) {
            first++;
            return logP;
        }

        if (first == 0 || !dynamics.areDynamicsKnown()) {
            precomputeMaxIntervalCache();
            dynamics.setDynamicsKnown();
        }

        if (!computeLikelihood) {
            first++;
            return 0;
        }

        // Running cursor over the precomputed cells. Advances monotonically
        // because currentTime advances monotonically through the tree.
        int rateCell = 0;

        do {
            double L = nextTreeEvent;

            if (L > 0) {
                int nSub = (Double.isFinite(maxInt) && L > maxInt)
                           ? (int) Math.ceil(L / maxInt) : 1;
                double subL = L / nSub;
                for (int k = 0; k < nSub; k++) {
                    double tMid = currentTime + (k + 0.5) * subL;
                    while (rateCell < mxNumCells - 1 && tMid > mxRateShiftEnds[rateCell]) {
                        rateCell++;
                    }
                    coalescentRates = mxCoalRates[rateCell];
                    double[] migrationRates = mxMigRates[rateCell];
                    logP += doEulerAtTime(subL, coalescentRates, migrationRates, mxIndicators);
                    if (logP == Double.NEGATIVE_INFINITY) {
                        first++;
                        return logP;
                    }
                }
                currentTime += L;
            }

            // Set rates at the event time so coalesce() uses the correct coalescentRates field.
            while (rateCell < mxNumCells - 1 && currentTime > mxRateShiftEnds[rateCell]) {
                rateCell++;
            }
            coalescentRates = mxCoalRates[rateCell];

            IntervalType type = treeIntervals.getIntervalType(treeInterval);
            if (type == IntervalType.COALESCENT) {
                nrLineages--;
                logP += coalesce(treeInterval, 0, 0.0, Double.POSITIVE_INFINITY);
            } else if (type == IntervalType.SAMPLE) {
                nrLineages++;
                sample(treeInterval, 0, 0.0, Double.POSITIVE_INFINITY);
            }

            if (logP == Double.NEGATIVE_INFINITY) {
                first++;
                return logP;
            }

            treeInterval++;
            try {
                nextTreeEvent = treeIntervals.getInterval(treeInterval);
            } catch (Exception e) {
                break;
            }
        } while (nextTreeEvent <= Double.POSITIVE_INFINITY);

        first++;
        return logP;
    }

    /**
     * Run a single Euler/RK4 step over `duration` with rates and indicators
     * supplied directly (no grid lookup). Mirrors doEuler() but bypasses the
     * ratesInterval index path.
     */
    public double doEulerAtTime(double duration, double[] coalRates, double[] migRates, int[] inds) {
        doEulerCallCount++;
        if (linProbs_tmp.length != linProbsLength + 1) {
            linProbs_tmp = new double[linProbsLength + 1];
        }
        System.arraycopy(linProbs, 0, linProbs_tmp, 0, linProbsLength);
        linProbs_tmp[linProbsLength] = 0;
        if (linProbsLength > 0) {
            linProbs[linProbsLength - 1] = 0;
        }

        // Use initWithIndicators when indicators are present so the integrator
        // can skip iterating over zero migration entries (sparsity win for
        // BSSVS-style models). The vendored Euler2ndOrder fork in
        // mascotdatastreams.ode patches the upstream bug where
        // initWithIndicators forgot to allocate sumDotStates. Fall back to
        // init when indicators are unavailable.
        if (inds != null && inds.length > 0) {
            euler.initWithIndicators(migRates, inds, coalRates, nrLineages);
        } else {
            euler.init(migRates, coalRates, nrLineages);
        }
        euler.calculateValues(duration, linProbs_tmp, linProbsLength + 1);

        System.arraycopy(linProbs_tmp, 0, linProbs, 0, linProbsLength);
        return linProbs_tmp[linProbsLength];
    }

	protected void setUpDynamics() {
    	int n = dynamics.getEpochCount();
    	double [][] coalescentRates = new double[n][];
    	double [][] migrationRates = new double[n][];
    	int [][] indicators = new int[n][];
    	double [] nextRateShift = dynamics.getIntervals();
    	for (int i = 0; i < n; i++) {
    		coalescentRates[i] = dynamics.getCoalescentRate(i);
            migrationRates[i] = dynamics.getBackwardsMigration(i);
    		indicators[i] = dynamics.getIndicators(i);
    	}
    	dynamics.setDynamicsKnown();
		euler.setUpDynamics(coalescentRates, migrationRates, indicators, nextRateShift);

	}

	double [] storedMigrationRates = new double[0];
    double [] storedCoalescentRates = new double[0];
    int storedNrLineages = -1;

	public double doEuler(double nextEventTime, int ratesInterval) {
		doEulerCallCount++;
		//for (int i = 0; i < linProbs.length; i++) linProbs_tmp[i] = linProbs[i];
		if (linProbs_tmp.length != linProbsLength + 1) {
			linProbs_tmp= new double[linProbsLength + 1];
		}
		System.arraycopy(linProbs,0,linProbs_tmp,0,linProbsLength);
		linProbs_tmp[linProbsLength] = 0;

		linProbs[linProbsLength-1] = 0;


//		if (dynamics.hasIndicators) {
//			euler.initWithIndicators(migrationRates, indicators, coalescentRates, nrLineages);
//			euler.calculateValues(nextEventTime, linProbs_tmp, linProbsLength + 1);
//		} else {
			euler.initAndcalculateValues(ratesInterval, nrLineages, nextEventTime, linProbs_tmp, linProbsLength + 1);
//		}

		//		System.out.println(Arrays.toString(linProbs));

		//for (int i = 0; i < linProbs.length; i++) linProbs[i] = linProbs_tmp[i];
		System.arraycopy(linProbs_tmp,0,linProbs,0,linProbsLength);
		return linProbs_tmp[linProbsLength];
	}


//    private void integrate(double duration){
//    }
//
//    private void ei(double duration, double[] linProbs_for_ode, double[] meanLinProbs){
//    	eulerIntegration(duration, linProbs_for_ode, meanLinProbs);
//    }
    
//    private double normalizeLineages(double [] linProbs){
//    	if (linProbs==null)
//    		return 0.0;
//
//
//    	double interval = 0.0;
//    	for (int i = 0; i < nrLineages; i++) {
//    		double lineProbs = 0.0;
//    		int u = i * states;
//    		for (int j = 0; j < states; j++) {
//    			if (linProbs[u]>=0.0){
//    				lineProbs += linProbs[u];
//    			} else {
//    				// try recalculation after lowering the tolerance
//    				System.out.println(linProbs[u]);
//    				recalculateLogP = true;
//    				return Math.log(1.0);
//    			}
//    			u++;
//    		}
//			if (lineProbs==0.0) {
//				return Double.NEGATIVE_INFINITY;
//			}
//			u = i * states;
//    		for (int j = 0; j < states; j++) {
//    			linProbs[u] = linProbs[u]/lineProbs;
//    			u++;
//    		}
//    		interval += lineProbs;
//    	}
//		// return mean P_t(T)
//		return Math.log(interval/(nrLineages));
//
//    }
//
    protected void sample(int currTreeInterval, int currRatesInterval, double nextTreeEvent, double nextRateShift) {
    	if (debug) {
    		System.err.println("sample activeLineages " + currTreeInterval + " = " + activeLineages);
    	}
		int incomingLines = treeIntervals.getLineagesAdded(currTreeInterval);
		int newLength = linProbsLength + 1 * states;

		int currPosition = linProbsLength;

		/*
		 * If there is no trait given as Input, the model will simply assume that
		 * the last value of the taxon name, the last value after a _, is an integer
		 * that gives the type of that taxon
		 */
		if (dynamics.typeTraitInput.get()!=null){
			Integer l = incomingLines; {
				activeLineages.add(l);//.getNr());
				int sampleState = nodeType[l];//dynamics.getValue(tree.getNode(l).getID());

				if (sampleState>= dynamics.getDimension()){
					System.err.println("sample discovered with higher state than dimension");
//					System.exit(1);
				}

				for (int i = 0; i < states; i++){
					if (i == sampleState){
						linProbs[currPosition] = 1.0;currPosition++;
					}
					else{
						linProbs[currPosition] = 0.0;currPosition++;
					}
				}
			}
		}else{
			Integer l = incomingLines; {
				activeLineages.add(l);//.getNr());
				String sampleID = tree.getNode(l).getID();
				int sampleState = 0;
				if (states > 1){
					String[] splits = sampleID.split("_");
					sampleState = Integer.parseInt(splits[splits.length-1]); //samples states (or priors) should eventually be specified in the XML
				}
				for (int i = 0; i < states; i++){
					if (i == sampleState){
						linProbs[currPosition] = 1.0;currPosition++;
					}
					else{
						linProbs[currPosition] = 0.0;currPosition++;
					}
				}
			}
		}
		linProbsLength = newLength;
		// store the node
       	storeNode(currTreeInterval, currRatesInterval, linProbs, logP, activeLineages, nextTreeEvent, nextRateShift, incomingLines);
    }

    protected double coalesce(int currTreeInterval, int currRatesInterval, double nextTreeEvent, double nextRateShift) {
    	int coalLines0 = treeIntervals.getLineagesRemoved(currTreeInterval,0);
    	int coalLines1 = treeIntervals.getLineagesRemoved(currTreeInterval,1);
    	if (debug) {
    		System.err.println("coalesce activeLineages " + currTreeInterval + " " + coalLines0 + " " + coalLines1 + " = " + activeLineages);
    	}

    	final int daughterIndex1 = activeLineages.indexOf(coalLines0);//.getNr());
		final int daughterIndex2 = activeLineages.indexOf(coalLines1);//.getNr());
		if (daughterIndex1 == -1 || daughterIndex2 == -1) {
			System.out.println(coalLines0/*.getNr()*/ + " " + coalLines1/*.getNr()*/ + " " + activeLineages);
			System.out.println("daughter lineages at coalescent event not found");
			System.exit(0);
			return Double.NaN;
		}
		double[] lambda = new double[states];
		double lambdaSum = 0;
		boolean isNegative = false;

		/*
		 * Calculate the overall probability for two strains to coalesce
		 * independent of the state at which this coalescent event is
		 * supposed to happen
		 */
		//double [] coalescentRates = ((Euler2ndOrder)euler).coalescentRates[Math.min(currRatesInterval, ((Euler2ndOrder)euler).coalescentRates[0].length - 1)];
        for (int k = 0; k < states; k++) {
        	Double pairCoalRate = coalescentRates[k] * linProbs[daughterIndex1*states + k] * linProbs[daughterIndex2*states + k];
			if (!Double.isNaN(pairCoalRate)){
				lambda[k] =  pairCoalRate;
				lambdaSum += pairCoalRate;
				if (pairCoalRate<0)
					isNegative = true;
			} else {
				return Double.NEGATIVE_INFINITY;
			}
        }

        int lineageToAdd = tree.getNode(coalLines0).getParent().getNr();
        activeLineages.add(lineageToAdd);

        // get the node state probabilities
		double[] pVec = new double[states];
		for (int i = 0; i < pVec.length; i++)
			pVec[i] = lambda[i]/lambdaSum;

		stateProbabilities[tree.getNode(coalLines0).getParent().getNr() - nrSamples] = pVec;

		//double[] linProbsNew  = new double[linProbsLength - states];

		int linCount = 0;
		// add all lineages execpt the daughter lineage to the new p array
		for (int i = 0; i <= nrLineages; i++){
			if (i != daughterIndex1 && i != daughterIndex2){
				for (int j = 0; j < states; j++){
					linProbsNew[linCount*states + j] = linProbs[i*states + j];
				}
				linCount++;
			}
		}
		// add the parent lineage
		for (int j = 0; j < states; j++){
			linProbsNew[linCount*states + j] = pVec[j];
		}
		// set p to pnew
		//double [] tmp = linProbs;
		linProbs = linProbsNew;
		linProbsNew = linProbs;
		linProbsLength = linProbsLength - states;


		//Remove daughter lineages from the line state probs
		if (daughterIndex1>daughterIndex2){
			// remove the daughter lineages from the active lineages
			activeLineages.remove(daughterIndex1);
			activeLineages.remove(daughterIndex2);
		} else {
			// remove the daughter lineages from the active lineages
			activeLineages.remove(daughterIndex2);
			activeLineages.remove(daughterIndex1);
		}

		if (isNegative){
			System.err.println("Coalescent probability is: " + Arrays.toString(lambda));
			return Double.NEGATIVE_INFINITY;
		}

		// store the node
        storeNode(currTreeInterval, currRatesInterval, linProbs, logP + Math.log(lambdaSum), activeLineages, nextTreeEvent, nextRateShift, lineageToAdd);

		if (lambdaSum==0)
			return Double.NEGATIVE_INFINITY;
		else
			return Math.log(lambdaSum);
    }


    public double[] getStateProb(int nr){
    	if (mascotImpl != null) {
    		return mascotImpl.getStateProb(nr);
    	}
    	return stateProbabilities[nr - nrSamples];
    }

    public double[] getRootState(){
    	if (mascotImpl != null) {
    		return mascotImpl.getRootState();
    	}
    	return stateProbabilities[stateProbabilities.length-1];
    }

    public String getType(){
   		return "state";
    }

    private void storeNode(int storingTreeInterval, int storingRatesInterval, double[] storeLinProbs,
		double probability, ArrayList<Integer> storeActiveLineages, double nextTreeEvent, double nextRateShift,
		int addedLineage) {
    	if (!useCache) {
    		return;
    	}

    	coalRatesInterval[storingTreeInterval] = storingRatesInterval;
    	int offset = 0;
    	if (storingTreeInterval > 0) {
    		offset = coalLinProbsLengths[storingTreeInterval-1];
    	}
    	System.arraycopy(storeLinProbs, 0, coalLinProbs, offset, linProbsLength);
    	coalLinProbsLengths[storingTreeInterval] = offset + linProbsLength;
    	coalLogP[storingTreeInterval] = probability;
    	nextTreeEvents[storingTreeInterval] = nextTreeEvent;
    	nextRateShifts[storingTreeInterval] = nextRateShift;
    }

    private int restoreNode(int restoringInterval){
    	//Log.warning("Restoring " + first + " " + restoringInterval);
    	int offset = 0;
    	if (restoringInterval > 0) {
    		offset = coalLinProbsLengths[restoringInterval-1];
    	}
    	linProbsLength = coalLinProbsLengths[restoringInterval] - offset;
    	System.arraycopy(coalLinProbs, offset, linProbs, 0, linProbsLength);

    	logP = coalLogP[restoringInterval];
    	return coalRatesInterval[restoringInterval + 1];
    }

    @Override
	public void store() {
    	if (!useCache) {
    		super.store();
    		return;
    	}
    	if (mascotImpl != null) {
    		mascotImpl.store();
    		super.store();
    		return;
    	}

    	storeLinP();
    	System.arraycopy(coalLogP, 0, storeLogP, 0, coalLogP.length);
    	System.arraycopy(coalRatesInterval, 0, storeRatesInterval, 0, coalRatesInterval.length);


    	System.arraycopy(nextTreeEvents, 0, storedNextTreeEvents, 0, nextTreeEvents.length);
    	System.arraycopy(nextRateShifts, 0, storedNextRateShifts, 0, nextRateShifts.length);


    	super.store();
    }

    private void storeLinP() {
    	System.arraycopy(coalLinProbsLengths, 0, storedCoalLinProbsLengths, 0, coalLinProbsLengths.length);
    	System.arraycopy(coalLinProbs, 0, storeLinProbs, 0, coalLinProbsLengths[coalLinProbsLengths.length - 1]);
//    	// store intermediate results
//    	for (int i = 0; i < coalLinProbs.length; i++) {
//    		double [] p = coalLinProbs[i];
//    		double [] q = storeLinProbs[i];
//    		if (p.length == q.length) {
//    			System.arraycopy(p, 0, q, 0, p.length);
//    		} else {
//    			q= Arrays.copyOf(p, p.length);
//    		}
//    	}
	}


	@Override
	public void restore(){
    	if (!useCache) {
    		super.restore();
    		return;
    	}

    	if (mascotImpl != null) {
    		mascotImpl.restore();
    		super.restore();
    		return;
    	}

    	// restore intermediate results
    	double [] tmp = storeLogP;
    	storeLogP = coalLogP;
    	coalLogP = tmp;

    	tmp = coalLinProbs;
    	coalLinProbs = storeLinProbs;
    	storeLinProbs = tmp;

    	int [] tmp2 = coalLinProbsLengths;
    	coalLinProbsLengths = storedCoalLinProbsLengths;
    	storedCoalLinProbsLengths = tmp2;

    	tmp2 = coalRatesInterval;
    	coalRatesInterval = storeRatesInterval;
    	storeRatesInterval = tmp2;

		tmp = nextTreeEvents;
		nextTreeEvents = storedNextTreeEvents;
		storedNextTreeEvents = tmp;

		tmp = nextRateShifts;
		nextRateShifts = storedNextRateShifts;
		storedNextRateShifts = tmp;

    	super.restore();
    }

    @Override
    protected boolean requiresRecalculation() {
        return ((CalculationNode) dynamics).isDirtyCalculation() || super.requiresRecalculation();
    }

}
