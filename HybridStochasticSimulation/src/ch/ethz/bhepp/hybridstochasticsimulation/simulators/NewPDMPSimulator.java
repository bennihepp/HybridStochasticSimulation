package ch.ethz.bhepp.hybridstochasticsimulation.simulators;
//package ch.ethz.khammash.hybridstochasticsimulation.simulators;
//
//import static com.google.common.base.Preconditions.checkArgument;
//
//import java.util.ArrayList;
//import java.util.Arrays;
//import java.util.HashSet;
//import java.util.LinkedList;
//import java.util.List;
//import java.util.Set;
//
//import org.apache.commons.math3.analysis.UnivariateFunction;
//import org.apache.commons.math3.analysis.solvers.BracketingNthOrderBrentSolver;
//import org.apache.commons.math3.analysis.solvers.UnivariateSolver;
//import org.apache.commons.math3.exception.MaxCountExceededException;
//import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
//import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;
//import org.apache.commons.math3.ode.sampling.FixedStepHandler;
//import org.apache.commons.math3.ode.sampling.StepHandler;
//import org.apache.commons.math3.ode.sampling.StepInterpolator;
//import org.apache.commons.math3.ode.sampling.StepNormalizer;
//import org.apache.commons.math3.ode.sampling.StepNormalizerBounds;
//import org.apache.commons.math3.random.RandomDataGenerator;
//import org.apache.commons.math3.util.FastMath;
//
//import ch.ethz.khammash.hybridstochasticsimulation.Utilities;
//import ch.ethz.khammash.hybridstochasticsimulation.math.MathUtilities;
//import ch.ethz.khammash.hybridstochasticsimulation.math.RandomDataUtilities;
//import ch.ethz.khammash.hybridstochasticsimulation.models.AdaptiveMSHRNModel;
//import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPMSHRNModel;
//import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
//import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
//import ch.ethz.khammash.hybridstochasticsimulation.networks.AdaptiveMSHRN;
//import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork;
//import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork.ReactionType;
//import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork.SpeciesType;
//import ch.ethz.khammash.hybridstochasticsimulation.simulators.ode.ExtendedStateObserverAdapter;
//import ch.ethz.khammash.hybridstochasticsimulation.simulators.ode.OdeAdapter;
//import ch.ethz.khammash.hybridstochasticsimulation.simulators.ode.StateObserverAdapter;
//import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;
//import ch.ethz.bhepp.ode.EventFunction;
//import ch.ethz.bhepp.ode.EventObserver;
//import ch.ethz.bhepp.ode.Solver;
//import ch.ethz.bhepp.ode.StateObserver;
//import ch.ethz.bhepp.ode.TimepointProvider;
//import ch.ethz.bhepp.ode.EventObserver.EventAction;
//import ch.ethz.bhepp.ode.lsodar.LsodarDirectSolver;
//
//public class NewPDMPSimulator extends AbstractSimulator<PDMPModel> {
//
//	public static final double DEFAULT_EVENT_HANDLER_MAX_CHECK_INTERVAL = Double.POSITIVE_INFINITY;
//	public static final double DEFAULT_EVENT_HANDLER_CONVERGENCE = 1e-12;
//	public static final double DEFAULT_EVENT_HANDLER_CONVERGENCE_FACTOR = Double.NaN;
//	public static final int DEFAULT_EVENT_HANDLER_MAX_ITERATION_COUNT = 1000;
//
//	private RandomDataGenerator rdg;
//	private List<TrajectoryRecorder> trajectoryRecorders;
//	private List<TrajectoryRecorder> optionalTrajectoryRecorders;
//	private List<TrajectoryRecorder> simulationInformationTrajectoryRecorders;
//	private Solver solver;
//
//	public NewPDMPSimulator() {
//		this(null, null);
//	}
//
//	public NewPDMPSimulator(RandomDataGenerator rdg) {
//		this(null, rdg);
//	}
//
//	public NewPDMPSimulator(Solver solver) {
//		this(solver, null);
//	}
//
//	public NewPDMPSimulator(Solver solver, RandomDataGenerator rdg) {
//		if (rdg == null)
//			rdg = new RandomDataGenerator();
//		this.rdg = rdg;
//		if (solver == null)
//			solver = LsodarDirectSolver.getInstance();
//		this.solver = solver;
//		trajectoryRecorders = new LinkedList<TrajectoryRecorder>();
//		optionalTrajectoryRecorders = new LinkedList<TrajectoryRecorder>();
//		simulationInformationTrajectoryRecorders = new LinkedList<TrajectoryRecorder>();
//	}
//
//	private class DefaultSimulationInformation implements SimulationInformation {
//    	private long integrationSteps = 0;
//    	private long reactionCount = 0;
//    	private boolean integrating = false;
//		@Override
//		public long getIntegrationSteps() {
//			return integrationSteps;
//		}
//		@Override
//		public long getReactionCount() {
//			return reactionCount;
//		}
//		@Override
//		public boolean isIntegrating() {
//			return integrating;
//		}
//		@Override
//		public double[] computeInformationState() {
//			double[] state = new double[3];
//			state[0] = integrating ? 1.0 : -1.0;
//			state[1] = integrationSteps;
//			state[2] = reactionCount;
//			return state;
//		}
//	}
//
//	private double stopTime;
//
////	private class Record {
////		public double t;
////		public double[] x;
////		public Record(double t, double[] x) {
////			this.t = t;
////			this.x = x.clone();
////		}
////	}
//
////	private List<Record> recordList;
//
//	final int NUM_OF_TIME_POINTS = 201;
////	int lastIndex;
//	double[] tSeries;
////	double[][] xSeries;
//	double[][] primaryStates;
////	ContinuousOutputModel com;
//
////	private class TemporaryDenseTrajectoryRecorder implements StepHandler {
////
////		private StepInterpolator[] interpolatorArray;
////		private double[] tSeries;
////		private int numOfTimepointsToUse;
////		private int lastIndexUsed;
////
////		public TemporaryDenseTrajectoryRecorder(int maxNumOfTimepoints) {
////			interpolatorArray = new StepInterpolator[maxNumOfTimepoints];
////			setNumOfTimepointsToUse(maxNumOfTimepoints);
////		}
////
////		public void setNumOfTimepointsToUse(int numOfTimepointsToUse) {
////			this.numOfTimepointsToUse = numOfTimepointsToUse;
////		}
////
////		public int findNearestLowerIndex(double t) {
////			int index = Arrays.binarySearch(tSeries, 0, lastIndexUsed, t);
////			if (index >= 0)
////				return index;
////			index = -index - 1;
////			return index;
////		}
////
////		public double[] getInterpolatedState(double t) {
////			int index = findNearestLowerIndex(t);
////			StepInterpolator interpolator = interpolatorArray[index];
////			interpolator.setInterpolatedTime(t);
////			return interpolator.getInterpolatedState();
////		}
////
////		public void rewind(double t) {
////			lastIndexUsed = findNearestLowerIndex(t);
////		}
////
////		@Override
////		public void init(double t0, double[] x0, double t1) {
////			this.lastIndexUsed = 0;
////		}
////
////		@Override
////		public void handleStep(StepInterpolator interpolator, boolean isLast) throws MaxCountExceededException {
////			if (lastIndexUsed < interpolatorArray.length) {
////				interpolatorArray[lastIndexUsed] = interpolator.copy();
////				tSeries[lastIndexUsed] = interpolator.getCurrentTime();
////				lastIndexUsed++;
////			} else
////				throw new RuntimeException("Too many calls to handleStep");
////		}
////
////	}
//
//	private class TemporaryTrajectoryRecorder implements StepHandler {
//
//		private double[][] xSeries;
//		private double t0;
//		private double t1;
//		private int numOfTimepointsToUse;
//		private int lastIndexUsed;
//		private double[] temp;
//		private MSHybridReactionNetwork network;
//		private double[] scaleFactors;
//
//		public TemporaryTrajectoryRecorder(int maxNumOfTimepoints, PDMPModel model) {
//			network = ((PDMPMSHRNModel)model).getNetwork();
//			xSeries = new double[maxNumOfTimepoints][];
//			for (int i=0; i < xSeries.length; i++)
//				xSeries[i] = new double[model.getVectorField().getDimension()];
//			temp = new double[model.getVectorField().getDimension()];
//			scaleFactors = new double[model.getVectorField().getDimension()];
//			setNumOfTimepointsToUse(maxNumOfTimepoints);
//		}
//
//		public void setNumOfTimepointsToUse(int numOfTimepointsToUse) {
//			this.numOfTimepointsToUse = numOfTimepointsToUse;
//		}
//
//		public int findNearestLowerIndex(double t) {
//			int index = (int)FastMath.floor((numOfTimepointsToUse - 1) * (t - t0) / (t1 - t0));
//			if (index < 0 || index >= numOfTimepointsToUse)
//				throw new IllegalArgumentException(String.format("Time t=%f is outside of the allowed range (%f,%f)", t, t0, t1));
//			return index;
//		}
//
//		public double[] getInterpolatedState(double t) {
//			double index = (numOfTimepointsToUse - 1) * (t - t0) / (t1 - t0);
//			int index1 = (int)FastMath.floor(index);
//			int index2 = (int)FastMath.ceil(index);
//			if (index1 == index2 || index2 >= numOfTimepointsToUse)
//				return xSeries[index1];
//			else {
//				double w1 = index - index1;
//				double w2 = index2 - index;
//				for (int s=0; s < xSeries[index1].length; s++) {
//					double x1 = xSeries[index1][s];
//					double x2 = xSeries[index2][s];
//					double v = x1 * w1 + x2 * w2;
//					temp[s] = v;
//				}
//				return temp;
//			}
//		}
//
//		public void rewind(double t) {
//			lastIndexUsed = findNearestLowerIndex(t) + 1;
//		}
//
//		@Override
//		public void init(double t0, double[] x0, double t1) {
//			this.t0 = t0;
//			this.t1 = t1;
//			this.lastIndexUsed = 0;
//			for (int s=0; s < scaleFactors.length - 1; s++)
//				scaleFactors[s] = network.getSpeciesScaleFactor(s);
//			scaleFactors[scaleFactors.length - 1] = 1.0;
////			handleStep(t0, x0, null, false);
//		}
//
////		@Override
////		public void handleStep(double t, double[] x, double[] xDot, boolean isLast) {
////			int index = findNearestLowerIndex(t);
////			while (lastIndexUsed < index) {
////				if (lastIndexUsed > 0) {
////					int prevIndex = lastIndexUsed - 1;
////					for (int s=0; s < x.length; s++)
////						xSeries[lastIndexUsed][s] = xSeries[prevIndex][s];
////				}
////				lastIndexUsed++;
////			}
////			if (lastIndexUsed == index) {
////				for (int s=0; s < x.length; s++)
////					xSeries[index][s] = x[s];
////				lastIndexUsed++;
////			}
////		}
//
//		@Override
//		public void handleStep(StepInterpolator interpolator, boolean isLast) throws MaxCountExceededException {
//			double t = interpolator.getCurrentTime();
//			int index = findNearestLowerIndex(t);
//			while (lastIndexUsed < index) {
//				double lastTime = t0 + (t1 - t0) * lastIndexUsed / (double)(numOfTimepointsToUse - 1); 
//				interpolator.setInterpolatedTime(lastTime);
//				double[] x = interpolator.getInterpolatedState();
//				for (int s=0; s < x.length; s++)
//					xSeries[lastIndexUsed][s] = x[s] * getScaleFactor(s);
//				lastIndexUsed++;
//			}
//			if (lastIndexUsed == index) {
//				double time = t0 + (t1 - t0) * index / (double)(numOfTimepointsToUse - 1);
//				interpolator.setInterpolatedTime(time);
//				double[] x = interpolator.getInterpolatedState();
//				for (int s=0; s < x.length; s++)
//					xSeries[index][s] = x[s] * getScaleFactor(s);
//				lastIndexUsed++;
//			}
//		}
//
//		private double getScaleFactor(int s) {
//			return scaleFactors[s];
//		}
//
//	}
//
//	private class RootFunction implements UnivariateFunction {
//
//		public UnivariateFunction function;
//    	public double offset;
//
//		@Override
//		public double value(double x) {
//			return function.value(x) - offset;
//		}
//
//	}
//
//	private class TemporaryFiniteTimepointProvider implements TimepointProvider {
//
//		private double t0;
//		private double t1;
//		private int numOfSteps;
//		private int index;
//
//		public void setLimits(double t0, double t1, int numOfSteps) {
//			this.t0 = t0;
//			this.t1 = t1;
//			this.numOfSteps = numOfSteps;
//			reset();
//		}
//
//		@Override
//		public double getInitialTimepoint() {
//			return t0;
//		}
//
//		@Override
//		public double getLastTimepoint() {
//			return t1;
//		}
//
//		@Override
//		public void reset() {
//			this.index = 0;
//		}
//
//		@Override
//		public double getCurrentTimepoint() {
//			double t = t0 + index / (double)(numOfSteps - 1) * (t1 - t0);
//			return t;
//		}
//
//		@Override
//		public boolean hasNextTimepoint(double tCurrent) {
//			return index + 1 < numOfSteps;
//		}
//
//		@Override
//		public double getNextTimepoint(double tCurrent) throws NoMoreTimepointsException {
//			index++;
//			return getCurrentTimepoint();
//		}
//		
//	}
//
//	private class TemporaryStateObserver implements StateObserver {
//
//		private double[][] xSeries;
//		private double t0;
//		private double t1;
//		private int numOfTimepointsToUse;
//		private int lastIndexUsed;
//		private double[] temp;
//		private MSHybridReactionNetwork network;
//		private double[] scaleFactors;
//
//		public TemporaryStateObserver(int maxNumOfTimepoints, PDMPMSHRNModel model) {
//			network = model.getNetwork();
//			xSeries = new double[maxNumOfTimepoints][];
//			for (int i=0; i < xSeries.length; i++)
//				xSeries[i] = new double[model.getVectorField().getDimension()];
//			temp = new double[model.getVectorField().getDimension()];
//			scaleFactors = new double[model.getVectorField().getDimension()];
//			setNumOfTimepointsToUse(maxNumOfTimepoints);
//		}
//
//		public void setNumOfTimepointsToUse(int numOfTimepointsToUse) {
//			this.numOfTimepointsToUse = numOfTimepointsToUse;
//		}
//
//		public int findNearestLowerIndex(double t) {
//			int index = (int)FastMath.floor((numOfTimepointsToUse - 1) * (t - t0) / (t1 - t0));
//			if (index < 0 || index >= numOfTimepointsToUse)
//				throw new IllegalArgumentException(String.format("Time t=%f is outside of the allowed range (%f,%f)", t, t0, t1));
//			return index;
//		}
//
//		public double[] getInterpolatedState(double t) {
//			double index = (numOfTimepointsToUse - 1) * (t - t0) / (t1 - t0);
//			int index1 = (int)FastMath.floor(index);
//			int index2 = (int)FastMath.ceil(index);
//			if (index1 == index2 || index2 >= numOfTimepointsToUse)
//				return xSeries[index1];
//			else {
//				double w1 = index - index1;
//				double w2 = index2 - index;
//				for (int s=0; s < xSeries[index1].length; s++) {
//					double x1 = xSeries[index1][s];
//					double x2 = xSeries[index2][s];
//					double v = x1 * w1 + x2 * w2;
//					temp[s] = v;
//				}
//				return temp;
//			}
//		}
//
//		private double getScaleFactor(int s) {
//			return scaleFactors[s];
//		}
//
//		@Override
//		public void initialize(double t0, double[] x0, double t1) {
//			this.t0 = t0;
//			this.t1 = t1;
//			this.lastIndexUsed = 0;
//			for (int s=0; s < scaleFactors.length - 1; s++)
//				scaleFactors[s] = network.getSpeciesScaleFactor(s);
//			scaleFactors[scaleFactors.length - 1] = 1.0;
////			handleStep(t0, x0, null, false);
//		}
//
//		@Override
//		public void report(double t, double[] x) {
//			int index = findNearestLowerIndex(t);
//			while (lastIndexUsed < index) {
//				if (lastIndexUsed > 0) {
//					int prevIndex = lastIndexUsed - 1;
//					for (int s=0; s < x.length; s++)
//						xSeries[lastIndexUsed][s] = xSeries[prevIndex][s];
//				}
//				lastIndexUsed++;
//			}
//			if (lastIndexUsed == index) {
//				for (int s=0; s < x.length; s++)
//					xSeries[index][s] = x[s] * getScaleFactor(s);
//				lastIndexUsed++;
//			}
//		}
//
//	}
//
//	public double simulate(final PDMPModel model, final double t0, final double[] x0, final double t1, double[] x1) {
////		rdg = new RandomDataGenerator();
////		rdg.reSeed(105L);
//
//		checkArgument(model.getNumberOfSpecies() == x0.length, "Expected model.getNumberOfSpecies() == x0.length but found %s != %s",
//				model.getNumberOfSpecies(), x0.length);
//
//		tSeries = MathUtilities.computeTimeSeries(NUM_OF_TIME_POINTS, t0, t1);
//
//    	// TODO
//    	// Compute dependency graph
//    	final AdaptiveMSHRNModel am = (AdaptiveMSHRNModel)model;
//    	final AdaptiveMSHRN net = (AdaptiveMSHRN)am.getNetwork();
//
//		double t = t0;
//		double[] x = new double[x0.length + 1];
//		for (int i=0; i < x0.length; i++)
//			x[i] = x0[i];
//		double[] xDot = new double[x.length];
//    	model.initialize(t, x);
//
////		xSeries = new double[x.length][];
////		for (int i=0; i < xSeries.length; i++) {
////			xSeries[i] = new double[NUM_OF_TIME_POINTS];
////		}
//		primaryStates = new double[NUM_OF_TIME_POINTS][];
//		for (int i=0; i < primaryStates.length; i++)
//			primaryStates[i] = new double[x0.length];
//
////    	// Simulation information
////		StateObserver stateObserver;
//		boolean recordSimulationInformation = simulationInformationTrajectoryRecorders.size() > 0 || printMessages;
//    	DefaultSimulationInformation simInfo = null;
//		if (recordSimulationInformation) {
//	    	simInfo = new DefaultSimulationInformation();
//	    	simInfo.integrationSteps = 0;
//	    	simInfo.reactionCount = 0;
//	    	simInfo.integrating = false;
////			stateObserver = new ExtendedStateObserverAdapter(simInfo, model, trajectoryRecorders, optionalTrajectoryRecorders, simulationInformationTrajectoryRecorders);
//		} else {
////			stateObserver = new StateObserverAdapter(model, trajectoryRecorders, optionalTrajectoryRecorders);
//		}
////		final StateObserver sto = stateObserver;
//////		recordList = new ArrayList<>(1000);
////		lastIndex = -1;
//		final TemporaryStateObserver tempStateObserver = new TemporaryStateObserver(1001, (PDMPMSHRNModel)model);
//		StateObserver stateObserver = new StateObserver() {
//
//			@Override
//			public void report(double t, double[] x) {
//				tempStateObserver.report(t, x);
//				int index = (int)FastMath.ceil((NUM_OF_TIME_POINTS - 1) * t / (t1 - t0));
//				record(index, t, x);
//			}
//
//			@Override
//			public void initialize(double t0, double[] x0, double t1) {
//				tempStateObserver.initialize(t0, x0, t1);
//			}
//
//			private void record(int index, double t, double[] x) {
//				double[] primaryState = model.computePrimaryState(t, x);
//				for (int s=0; s < primaryState.length; s++)
//					primaryStates[index][s] = primaryState[s];
//			}
//
//		};
////		stateObserver.initialize(t, x, t1);
////		List<PDMPEventObserver> optionalEventObservers = model.getOptionalEventObservers();
////		List<PDMPEventObserver> eventObservers = new ArrayList<PDMPEventObserver>(1 + optionalEventObservers.size());
////		eventObservers.add(model.getJumpEventObserver());
////		eventObservers.addAll(optionalEventObservers);
////		OdeAdapter ode = new OdeAdapter(model);
//////		EventObserverAdapter eventObserver = new EventObserverAdapter(eventObservers);
////		EventObserverAdapter eventObserver = new EventObserverAdapter(eventObservers) {
////
////			@Override
////			public void computeEventValues(double t, double[] x, double[] values) {
////				for (int i=0; i < values.length; i++)
////					values[i] = 1.0;
////			}
////
////		};
////		eventObserver.initialize(t, x, t1);
//		StochasticReactionNetworkModel rnm = model.getTransitionMeasure();
//    	FirstOrderDifferentialEquations vectorField = model.getVectorField();
//
//    	TemporaryFiniteTimepointProvider timepointProvider = new TemporaryFiniteTimepointProvider();
//
////    	FiniteTimepointProvider timepointProvider = null;
////		if (trajectoryRecorders.size() > 0) {
////			for (TrajectoryRecorder tr : trajectoryRecorders)
////				if (tr instanceof FiniteTrajectoryRecorder) {
////					FiniteTrajectoryRecorder ftr = (FiniteTrajectoryRecorder)tr; 
//////					timepointProvider = new FiniteTimepointProvider(ftr.gettSeries());
////					timepointProvider = new FiniteTimepointProvider(ftr.gettSeries()) {
////
////						@Override
////						public boolean hasNextTimepoint(double tCurrent) {
////							if (tCurrent >= stopTime)
////								return false;
////							return super.hasNextTimepoint(tCurrent);
////						}
////
////						@Override
////						public double getNextTimepoint(double tCurrent) throws NoMoreTimepointsException {
////							double nextTimepoint = super.getNextTimepoint(tCurrent);
////							if (nextTimepoint > stopTime)
////								return stopTime;
////							return nextTimepoint;
////						}
////
////					};
////					break;
////				}
////		}
////		if (timepointProvider == null)
////			timepointProvider = new FiniteTimepointProvider(t0, t1);
//
//		double[] propVec = new double[rnm.getNumberOfReactions()];
//    	double[] reactionCounterArray = new double[model.getNumberOfReactions()];
//    	double msgDt = (t1 - t0) / 20.0;
//    	double nextMsgT = t0 + msgDt;
//    	int j = 0;
//
//        double[] q = new double[x0.length];
//        double[] p = new double[x0.length];
//
//
//    	final UnivariateFunction propensityFunction = new UnivariateFunction() {
//
//			@Override
//			public double value(double t) {
////				double[] state = ttr.getInterpolatedState(t);
//				double[] state = tempStateObserver.getInterpolatedState(t);
//				return state[state.length - 1];
////				com.setInterpolatedTime(t);
////				double[] state = com.getInterpolatedState();
////				return state[state.length - 1];
////				if (t == tt)
////					return 0.0;
////				int index = (int)FastMath.floor((NUM_OF_TIME_POINTS - 1) * t / (t1 - t0));
////				if (index + 1 >= tSeries.length)
////					return xSeries[xSeries.length - 1][tSeries.length - 1];
////				if (index + 1 > lastIndex)
////					return xSeries[xSeries.length - 1][lastIndex];
////				double x1 = xSeries[xSeries.length - 1][index];
////				double x2 = xSeries[xSeries.length - 1][index + 1];
////				double t1 = tSeries[index];
////				double t2 = tSeries[index + 1];
////				double w1 = (t - t1) / (t2 - t1);
////				double w2 = 1.0 - w1;
////				double v = w1 * x1 + w2 * x2;
////				return v;
//			}
//
//		};
//
//		RootFunction rootFunction = new RootFunction();
//		rootFunction.function = propensityFunction;
//    	UnivariateSolver rootSolver = new BracketingNthOrderBrentSolver();
//
//    	Set<Integer> coupledDiscreteSpecies = new HashSet<>(net.getNumberOfSpecies());
//    	Set<Integer> coupledStochasticReactions = new HashSet<>(net.getNumberOfReactions());
//
//		OdeAdapter ode = new OdeAdapter(model);
//    	EventObserver eventObserver = new EventObserver() {
//
//			@Override
//			public EventAction report(int eventIndex, double t, double[] x) {
//				return EventAction.CONTINUE;
//			}
//
//			@Override
//			public void initialize(double t0, double[] x0, double t1) {
//			}
//
//		};
//		EventFunction eventFunction = new EventFunction() {
//
//			@Override
//			public int getNumberOfEventValues() {
//				return 0;
//			}
//
//			@Override
//			public void computeEventValues(double t, double[] x, double[] values) {
//			}
//
//		};
//
//		synchronized (solver) {
//			solver.initialize(ode, eventFunction, stateObserver, eventObserver);
//
//			final long startTime = System.currentTimeMillis();
//			while (true) {
//	
//				double tPrevious = t;
//	//			double[] xPrevious = x.clone();
//	
//				if (showProgress)
//					while (t > nextMsgT) {
//						System.out.println("Progress: " + (100 * (nextMsgT - t0)/(t1 - t0)) + "%");
//						nextMsgT += msgDt;
//					}
//				boolean hasDeterministicPart = model.hasVectorField();
//				// TODO: Check whether this gives any performance gain
//				if (hasDeterministicPart && model.isTimeIndependent()) {
//					vectorField.computeDerivatives(t, x, xDot);
//					boolean allZero = true;
//					for (int s=0; s < x0.length; s++)
//						if (xDot[s] != 0.0) {
//							allZero = false;
//							break;
//						}
//						if (allZero)
//							hasDeterministicPart = false;
//				}
//				double propSum = 0.0;
//				boolean propensitiesComputed = false;
//				if (hasDeterministicPart) {
//	
//		        	// TODO
//					coupledDiscreteSpecies.clear();
//		        	for (int r=0; r < net.getNumberOfReactions(); r++) {
//						if (net.getReactionType(r) == ReactionType.DETERMINISTIC) {
//							List<Integer> coupledSpeciesList = getSpeciesInfluencingReaction(net, r);
//							for (int s : coupledSpeciesList) {
//								if (net.getSpeciesType(s) == SpeciesType.DISCRETE)
//									coupledDiscreteSpecies.add(s);
//							}
//						}
//		        	}
//		        	coupledStochasticReactions.clear();
//		        	for (int s : coupledDiscreteSpecies) {
//		        		List<Integer> coupledReactionList = getReactionsChangingSpecies(net, s);
//		        		coupledStochasticReactions.addAll(coupledReactionList);
//		        	}
//		        	//
//		        	double coupledPropensitiesSum = 0.0;
//		        	for (int r : coupledStochasticReactions)
//		        		coupledPropensitiesSum += am.computePropensity(r, t, x);
//		        	double integrationTimeStep;
//		        	final double alpha = 0.90;
//		        	if (coupledPropensitiesSum > 0.0)
//			        	integrationTimeStep = - FastMath.log(1 - alpha) / coupledPropensitiesSum;
//	//	        		integrationTimeStep = 1 / coupledPropensitiesSum;
//		        	else
//		        		integrationTimeStep = 1 * (t1 - t);
//		        	//
//		        	stopTime = FastMath.min(t + integrationTimeStep, t1);
//		        	//
//	
//				}
//	
//				if (hasDeterministicPart) {
//	
//		        	if (recordSimulationInformation)
//		        		simInfo.integrating = true;
//	
//					// Evolve ODE until next stochastic reaction fires
//	//	        	double nextUnitJumpTime = -FastMath.log(rdg.nextUniform(0.0,  1.0));
//	//	        	// x[x.length - 1] will be < 0
//	//		        x[x.length - 1] = -nextUnitJumpTime;
//			        x[x.length - 1] = 0.0;
//			        model.checkAndHandleOptionalEvent(t, x);
//	//			        System.out.println(String.format("Evolving from t=%f to t=%f", t, stopTime));
////		        	ttr.setNumOfTimepointsToUse(101);
//	//		        double[] pm1 = model.computePrimaryState(t, x).clone();
////			        t = integrator.integrate(vectorField, t, x, stopTime, x);
//	//		        double[] pm2 = model.computePrimaryState(t, x).clone();
//					stateObserver.initialize(t, x, stopTime);
//			        timepointProvider.setLimits(t, stopTime, 101);
//		        	solver.prepare(timepointProvider, x);
//	        		try {
//	        			t = solver.integrate();
//	        		} catch (Exception e) {
//	        			throw e;
//	        		}
////		        	while (true) {
////		        		try {
////		        			t = solver.integrate();
////		        		} catch (Exception e) {
////		        			throw e;
////		        		}
////			        	if (recordSimulationInformation)
////			        		simInfo.integrationCount++;
////			        	if (showProgress)
////							while (t > nextMsgT) {
////								System.out.println("Progress: " + (100 * (nextMsgT - t0)/(t1 - t0)) + "%");
////								nextMsgT += msgDt;
////							}
////			        	if (model.hasOptionalEventOccured()) {
////				        	model.handleOptionalEvent(t, x);
////				        	timepointProvider.setCurrentTimepoint(t);
////				        	solver.prepare(timepointProvider, x);
////			        	} else
////			        		break;
////			        }
//		        	if (recordSimulationInformation)
//		        		simInfo.integrating = false;
//	
//		        	// TODO
//	//		        System.out.println(String.format("Evolved from t=%f to t=%f", tPrevious, t));
//			        // Find all reactions that fired between tPrevious and t
//					double nextUnitJumpTime = 0.0;
//			        double tReaction = tPrevious;
//			        final double tt = tPrevious;
//			        Arrays.fill(q, 0.0);
//			        int coupledReactionCounter = 0;
//			        int unCoupledReactionCounter = 0;
//					while (tReaction < t) {
//			        	nextUnitJumpTime += -FastMath.log(rdg.nextUniform(0.0,  1.0));
//			        	rootFunction.offset = nextUnitJumpTime;
//	//			        	final UnivariateInterpolator interpolator = new LinearInterpolator();
//	//			        	final UnivariateFunction propensityFunction = interpolator.interpolate(tSeries, xSeries[xSeries.length - 1]);
//	//			        	System.out.println(String.format("  rootFunction(tReaction)=%f, propensityFunction(tReaction)=%f, nextUnitJumpTime=%f", rootFunction.value(tReaction), propensityFunction.value(tReaction), nextUnitJumpTime));
//	//			        	System.out.println(String.format("  rootFunction(t)=%f, propensityFunction(t)=%f, nextUnitJumpTime=%f", rootFunction.value(t), propensityFunction.value(t), nextUnitJumpTime));
//	//			        	UnivariateFunction propensityFunction = new UnivariateFunction() {
//	//
//	//			        		private int findIndex(double t) {
//	//			        			int index = Collections.binarySearch(recordList, new Record(t, new double[0]), new Comparator<Record>() {
//	//
//	//									@Override
//	//									public int compare(Record o1, Record o2) {
//	//										return Double.compare(o1.t, o2.t);
//	//									}
//	//
//	//			        			});
//	//			        			if (index >= 0)
//	//			        				return index;
//	//			        			else
//	//			        				return -index - 1;
//	//			        		}
//	//
//	//							@Override
//	//							public double value(double t) {
//	//								int index = findIndex(t);
//	//								if (index + 1 < recordList.size()) {
//	//									double[] x1 = recordList.get(index).x;
//	//									double[] x2 = recordList.get(index + 1).x;
//	//									double t1 = recordList.get(index).t;
//	//									double t2 = recordList.get(index + 1).t;
//	//									double w1 = (t - t1) / (t2 - t1);
//	//									double w2 = 1.0 - t1;
//	//									double v = w1 * x1[x1.length - 1] + w2 * x2[x2.length - 2] - nextUnitJumpTime;
//	//									index++;
//	//									return v - nextUnitJumpTime;
//	//								} else if (index < recordList.size()) {
//	//									double[] x1 = recordList.get(index).x;
//	//									double v = x1[x1.length - 1];
//	//									return v - nextUnitJumpTime;
//	//								} else
//	//									return 0.0;
//	//							}
//	//		
//	//						};
//			        	try {
//							if (rootFunction.value(t) < 0) {
//	//								System.out.println(String.format("  No stochastic reaction took place: rootFunction(t)=%f", rootFunction.value(t)));
//								break;
//							}
//			        	} catch (Exception e) {
//			        		rootFunction.value(t);
//			        	}
//	
//			        	try {
//			        		tReaction = rootSolver.solve(1000, rootFunction, tReaction, t);
//	//			        		System.out.println(String.format("  found root at tReaction=%f, rootFunction(tReaction)=%f, propensityFunction(tReaction)=%f", tReaction, rootFunction.value(tReaction), propensityFunction.value(tReaction)));
//			        		nextUnitJumpTime = propensityFunction.value(tReaction);
//			        	} catch (Exception e) {
//	//				        	System.out.println(String.format("  rootFunction(tReaction)=%f, propensityFunction(tReaction)=%f, nextUnitJumpTime=%f", rootFunction.value(tReaction), propensityFunction.value(tReaction), nextUnitJumpTime));
//	//				        	System.out.println(String.format("  rootFunction(t)=%f, propensityFunction(t)=%f, nextUnitJumpTime=%f", rootFunction.value(t), propensityFunction.value(t), nextUnitJumpTime));
//			        		tReaction = rootSolver.solve(1000, rootFunction, tReaction, t);
//			        	}
////			        	double[] tmp = ttr.getInterpolatedState(tReaction);
//			        	double[] tmp = tempStateObserver.getInterpolatedState(tReaction);
//	//		        	computeInterpolatedScaledState(tReaction, x, net);
//			        	for (int i=0; i < q.length; i++)
//			        		p[i] = (tmp[i] + q[i]) * net.getInverseSpeciesScaleFactor(i);
//				        rnm.computePropensities(tReaction, p, propVec);
//				        int reaction = RandomDataUtilities.sampleFromProbabilityMassFunction(rdg, propVec);
//	//						System.out.println(String.format("  Stochastic reaction r=%d took place at t=%f", reaction, tReaction));
//				        if (coupledStochasticReactions.contains(reaction)) {
//	//				        	System.out.println("  Stochastic reaction was coupled to vector field");
//				        	t = tReaction;
//	//						System.out.println(String.format("  Rewinding to t=%f", t));
//				        	coupledReactionCounter++;
//				        } else {
//				        	unCoupledReactionCounter++;
//				        }
//				        rnm.changeState(reaction, t, q);
//			        	if (recordSimulationInformation)
//			        		simInfo.reactionCount++;
//			        	reactionCounterArray[reaction]++;
//	//			        	stateObserver.report(t, x);
//			        	// TODO: Make this value configurable (coupled to N?)
//					}
//	
//	//				System.out.println(String.format("  total=%d, coupled=%d, uncoupled=%d", coupledReactionCounter + unCoupledReactionCounter, coupledReactionCounter, unCoupledReactionCounter));
//	
//	//				computeInterpolatedScaledState(t, x, net);
//	//	        	for (int i=0; i < q.length; i++)
//	//	        		x[i] = tmp[i] + q[i];
//	
//					if (coupledReactionCounter > 0) {
////						double[] tmp = ttr.getInterpolatedState(t);
//						double[] tmp = tempStateObserver.getInterpolatedState(t);
//			        	for (int i=0; i < q.length; i++) {
//			        		x[i] = tmp[i] + q[i];
//			        		x[i] *= net.getInverseSpeciesScaleFactor(i);
//			        	}
//					} else {
//			        	for (int i=0; i < q.length; i++) {
//			        		x[i] += q[i] * net.getInverseSpeciesScaleFactor(i);
//			        	}
//					}
//	
////		        	stepHandler.handleStep(t, x, null, false);
//	
//				} else {
//	
//			        rnm.computePropensities(t, x, propVec);
//			        for (int i=0; i < propVec.length; i++)
//			        	propSum += propVec[i];
//			        propensitiesComputed = true;
//			        // Find next reaction time point
//			        if (propSum <= 0)
//			        	break;
//			        // -FastMath.log(rng.nextUniform(0.0,  1.0))
//			        double tau = rdg.nextExponential(1 / propSum);
//			        t = t + tau;
//	
//			        // Stop if we reached the end-timepoint
//			        if (t >= t1)
//			        	break;
//	
//			        // Determine which reaction fired and update state
//			        if (!propensitiesComputed) {
//				        rnm.computePropensities(t, x, propVec);
//				        for (int i=0; i < propVec.length; i++)
//				        	propSum += propVec[i];
//			        }
//			        int reaction = RandomDataUtilities.sampleFromProbabilityMassFunction(rdg, propVec);
//		    		rnm.changeState(reaction, t, x);
//	//			        double u = rdg.nextUniform(0.0, 1.0);
//	//			        double w = 0.0;
//	//			        int reaction = -1;
//	//			        for (int l=0; l < propVec.length; l++) {
//	//			        	w = w + propVec[l] / propSum;
//	//			        	if (u < w) {
//	//			        		reaction = l;
//	//			        		rnm.changeState(reaction, t, x);
//	//			        		break;
//	//			        	}
//	//			        }
//			        if (reaction >= 0) {
//			        	if (recordSimulationInformation)
//			        		simInfo.reactionCount++;
//			        	reactionCounterArray[reaction]++;
//	//			        	stateObserver.report(t, x);
////			        	stepHandler.handleStep(t, x, null, false);
//			        	// TODO: Make this value configurable (coupled to N?)
//			        	if (j > 100) {
//			        		model.checkAndHandleOptionalEvent(t, x);
//			        		j = 0;
//			        	}
//			        	j++;
//			        }
//				}
//	
//		        // Stop if we reached the end-timepoint
//		        if (t >= t1)
//		        	break;
//	
//			}
//
//			if (printMessages) {
//				System.out.println("Integrator invocations: " + simInfo.getIntegrationSteps());
//	//				System.out.println("Observed " + eventObserver.getEventCount() + " events");
//				final long endTime = System.currentTimeMillis();
//				System.out.println("Execution time: " + (endTime - startTime));
//		    	long evaluationCounter = ode.getEvaluations();
////				long evaluationCounter = integrator.getEvaluations();
//				System.out.println("Total of " + evaluationCounter + " evaluations and " + simInfo.getReactionCount() + " reactions performed");
//				Utilities.printArray("Total reaction counts", reactionCounterArray);
//				for (int r=0; r < reactionCounterArray.length; r++)
//					reactionCounterArray[r] /= simInfo.getReactionCount();
//				Utilities.printArray("Relative reaction counts", reactionCounterArray);
//			}
//			for (int i=0; i < x1.length; i++)
//				x1[i] = x[i];
//	//	    	for (TrajectoryRecorder<T> handler : trajectoryRecorders)
//	//	    		handler.setFinalState(t1, x1);
//	//        	stateObserver.report(t1, x1);
////	    	stepHandler.handleStep(t1, x1, null, false);
//	    	solver.dispose();
//	//	    	sto.initialize(t0, x0, t1);
//	//	    	for (int i=0; i < tSeries.length; i++) {
//	//	    		sto.report(t, primaryStates[i]);
//	//	    		double[] xx = new double[xSeries.length];
//	//	    		for (int k=0; k < xSeries.length; k++)
//	//	    			xx[k] = xSeries[k][i];
//	//	    		sto.report(tSeries[i], xx.clone());
//	//	    	}
//	    	for (TrajectoryRecorder tr : trajectoryRecorders) {
//	    		tr.beginRecording(t0, x0, t1);
//	    		for (int i=0; i < tSeries.length; i++)
//	    			tr.record(tSeries[i], primaryStates[i]);
//	    		tr.endRecording(primaryStates[primaryStates.length - 1]);
//	    	}
//		}
//
//		return t;
//	}
//
////	// Compute the interpolated scaled state at time t from the tSeries and primaryStates arrays
////	private void computeInterpolatedScaledState(double t, double[] x, MSHybridReactionNetwork net) {
//////		computeInterpolatedPrimaryState(t, x, net);
////		com.setInterpolatedTime(t);
////		double[] state = com.getInterpolatedState();
////		for (int s=0; s < net.getNumberOfSpecies(); s++)
////			x[s] = state[s] * net.getInverseSpeciesScaleFactor(s);
////	}
//
////	// Compute the interpolated primary state at time t from the tSeries and primaryStates arrays (scaled)
////	private void computeInterpolatedPrimaryState(double t, double[] x, MSHybridReactionNetwork net) {
////		double t0 = tSeries[0];
////		double t1 = tSeries[tSeries.length - 1];
////		int index = (int)FastMath.floor((tSeries.length - 1) * t / (t1 - t0));
////		if (index + 1 <= lastIndex) {
////			double tempT1 = tSeries[index];
////			double tempT2 = tSeries[index + 1];
////			double w1 = (t - tempT1) / (tempT2 - tempT1);
////			double w2 = 1.0 - w1;
////	    	for (int s=0; s < net.getNumberOfSpecies(); s++) {
////				double tempX1 = primaryStates[index][s];
////				double tempX2 = primaryStates[index + 1][s];
////				double v;
////				if (net.getSpeciesType(s) == SpeciesType.DISCRETE)
////					v = tempX2;
////				else
////					v = w1 * tempX1 + w2 * tempX2;
////				x[s] = v;
////	    	}
////		} else {
////	    	for (int s=0; s < net.getNumberOfSpecies(); s++) {
////	    		x[s] = primaryStates[lastIndex][s];
////	    	}
////		}
////	}
//
//	private List<Integer> getSpeciesInfluencingReaction(AdaptiveMSHRN net, int r) {
//		List<Integer> influcencingSpecies = new ArrayList<>(net.getNumberOfSpecies());
//		int[] consumptionStoichiometries = net.getConsumptionStoichiometries(r);
//		for (int s=0; s < consumptionStoichiometries.length; s++)
//			if (consumptionStoichiometries[s] > 0)
//				influcencingSpecies.add(s);
//		return influcencingSpecies;
//	}
//
//	private List<Integer> getReactionsChangingSpecies(MSHybridReactionNetwork net, int species) {
//		List<Integer> changingReactions = new ArrayList<>(net.getNumberOfReactions());
//		for (int r=0; r < net.getNumberOfReactions(); r++) {
//			if (net.getStoichiometry(species, r) != 0)
//				changingReactions.add(r);
//		}
//		return changingReactions;
//	}
//
//	@Override
//	public void addTrajectoryRecorder(TrajectoryRecorder tr) {
//		trajectoryRecorders.add(tr);
//	}
//
//	@Override
//	public void removeTrajectoryRecorder(TrajectoryRecorder tr) {
//		trajectoryRecorders.remove(tr);
//	}
//
//	@Override
//	public void clearTrajectoryRecorders() {
//		trajectoryRecorders.clear();
//	}
//
//	public void addOptionalTrajectoryRecorder(TrajectoryRecorder tr) {
//		optionalTrajectoryRecorders.add(tr);
//	}
//
//	public void removeOptionalTrajectoryRecorder(TrajectoryRecorder tr) {
//		optionalTrajectoryRecorders.remove(tr);
//	}
//
//	public void clearOptionalTrajectoryRecorders() {
//		optionalTrajectoryRecorders.clear();
//	}
//
//	public void addSimulationInformationTrajectoryRecorder(TrajectoryRecorder tr) {
//		simulationInformationTrajectoryRecorders.add(tr);
//	}
//
//	public void removeSimulationInformationTrajectoryRecorder(TrajectoryRecorder tr) {
//		simulationInformationTrajectoryRecorders.remove(tr);
//	}
//
//	public void clearSimulationInformationTrajectoryRecorders() {
//		simulationInformationTrajectoryRecorders.clear();
//	}
//
//}
