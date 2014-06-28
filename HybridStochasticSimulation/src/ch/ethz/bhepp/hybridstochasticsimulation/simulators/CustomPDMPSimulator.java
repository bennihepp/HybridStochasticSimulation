package ch.ethz.bhepp.hybridstochasticsimulation.simulators;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.BracketingNthOrderBrentSolver;
import org.apache.commons.math3.analysis.solvers.UnivariateSolver;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.util.FastMath;

import ch.ethz.bhepp.hybridstochasticsimulation.Utilities;
import ch.ethz.bhepp.hybridstochasticsimulation.math.distributions.DiscreteProbabilityDistribution;
import ch.ethz.bhepp.hybridstochasticsimulation.models.AdaptiveMSHRNModel;
import ch.ethz.bhepp.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.bhepp.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.AdaptiveMSHRN;
import ch.ethz.bhepp.hybridstochasticsimulation.simulators.ode.EventObserverAdapter;
import ch.ethz.bhepp.hybridstochasticsimulation.simulators.ode.ExtendedStateObserverAdapter;
import ch.ethz.bhepp.hybridstochasticsimulation.simulators.ode.OdeAdapter;
import ch.ethz.bhepp.hybridstochasticsimulation.simulators.ode.StateObserverAdapter;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.TrajectoryRecorder;
import ch.ethz.bhepp.ode.Solver;
import ch.ethz.bhepp.ode.lsodar.LsodarDirectSolver;
import ch.ethz.bhepp.ode.nonstiff.EulerSolver;

public class CustomPDMPSimulator extends AbstractSimulator<PDMPModel> {

	public static final double DEFAULT_EVENT_HANDLER_MAX_CHECK_INTERVAL = Double.POSITIVE_INFINITY;
	public static final double DEFAULT_EVENT_HANDLER_CONVERGENCE = 1e-12;
	public static final double DEFAULT_EVENT_HANDLER_CONVERGENCE_FACTOR = Double.NaN;
	public static final int DEFAULT_EVENT_HANDLER_MAX_ITERATION_COUNT = 1000;

	private RandomDataGenerator rdg;
	private List<TrajectoryRecorder> trajectoryRecorders;
	private List<TrajectoryRecorder> optionalTrajectoryRecorders;
	private List<TrajectoryRecorder> simulationInformationTrajectoryRecorders;
	private Solver solver;

	public CustomPDMPSimulator() {
		this(null, null);
	}

	public CustomPDMPSimulator(RandomDataGenerator rdg) {
		this(null, rdg);
	}

	public CustomPDMPSimulator(Solver solver) {
		this(solver, null);
	}

	public CustomPDMPSimulator(Solver solver, RandomDataGenerator rdg) {
		if (rdg == null)
			rdg = new RandomDataGenerator();
		this.rdg = rdg;
		if (solver == null)
			solver = LsodarDirectSolver.getInstance();
		this.solver = solver;
		trajectoryRecorders = new LinkedList<TrajectoryRecorder>();
		optionalTrajectoryRecorders = new LinkedList<TrajectoryRecorder>();
		simulationInformationTrajectoryRecorders = new LinkedList<TrajectoryRecorder>();
	}

	private class TemporaryTrajectoryInterpolator {

		private int maxNumOfTimepoints;
		private int dimension;
		private double[] tSeries;
		private double[][] xSeries;
		private double[] xTemp;
		private int nextIndex;
//		private double[][] qSeries;

		public TemporaryTrajectoryInterpolator(int maxNumOfTimepoints, int dimension) {
			this.maxNumOfTimepoints = maxNumOfTimepoints;
			this.dimension = dimension;
			tSeries = new double[maxNumOfTimepoints];
			xSeries = new double[maxNumOfTimepoints][];
			for (int i=0; i < maxNumOfTimepoints; i++)
				xSeries[i] = new double[dimension];
			xTemp = new double[dimension];
//			qSeries = new double[dimension][];
//			for (int i=0; i < dimension; i++)
//				qSeries[i] = new double[maxNumOfTimepoints];
		}

		public int getNumberOfRecordedTimepoints() {
			return nextIndex;
		}

		public double getRecordedTimepoint(int index) {
			return tSeries[index];
		}

		public double[] getRecordedState(int index) {
			return xSeries[index];
		}

		public int binarySearch(double t) {
			int index = Arrays.binarySearch(tSeries, 0, nextIndex, t);
			return index;
		}

		public double computeLeftWeight(double t, int index1, int index2) {
			double w = (tSeries[index2] - t) / (tSeries[index2] - tSeries[index1]);
			if (w < 0.0 || w > 1.0)
				throw new IllegalArgumentException();
			return w;
		}

		public double[] getInterpolatedState(double t) {
			int index = binarySearch(t);
			if (index >= 0) {
				return xSeries[index];
			} else {
				int index1 = (-index - 1) - 1;
				int index2 = index1 + 1;
				double w1 = computeLeftWeight(t, index1, index2);
				double w2 = 1.0 - w1;
				for (int i=0; i < dimension; i++) {
					double x1 = xSeries[index1][i];
					double x2 = xSeries[index2][i];
					double v = x1 * w1 + x2 * w2;
					xTemp[i] = v;
				}
				return xTemp;
			}
		}

		public double getInterpolatedState(double t, int i) {
			int index = binarySearch(t);
			if (index >= 0) {
				return xSeries[index][i];
			} else {
				int index1 = (-index - 1) - 1;
				int index2 = index1 + 1;
				double w1 = computeLeftWeight(t, index1, index2);
				double w2 = 1.0 - w1;
				double x1 = xSeries[index1][i];
				double x2 = xSeries[index2][i];
				double t1 = tSeries[index1];
				double t2 = tSeries[index2];
				double v = x1 * w1 + x2 * w2;
				return v;
			}
		}

		public void reset(double t0, double[] x0) {
			nextIndex = 0;
			record(t0, x0);
		}

		private void record(double t, double[] x) {
			if (nextIndex >= maxNumOfTimepoints)
				throw new ArrayIndexOutOfBoundsException();
			tSeries[nextIndex] = t;
			for (int i=0; i < dimension; i++)
				xSeries[nextIndex][i] = x[i];
//			for (int i=0; i < dimension; i++)
//				qSeries[i][nextIndex] = x[i];
			nextIndex++;
		}

	}

	private class InnerTrajectoryRecorder {

		private int numOfTimepoints;
		private double t0;
		private double t1;
		private int dimension;
		public boolean[] mask;
		public double[] tSeries;
		public double[][] xSeries;
//		private int nextIndex;

		public InnerTrajectoryRecorder(int numOfTimepoints) {
			this.numOfTimepoints = numOfTimepoints;
		}

		public void init(double t0, double[] x0, double t1) {
//			nextIndex = 0;
			this.t0 = t0;
			this.t1 = t1;
			dimension = x0.length;
			mask = new boolean[numOfTimepoints];
			tSeries = new double[numOfTimepoints];
			xSeries = new double[numOfTimepoints][];
			for (int i=0; i < numOfTimepoints; i++)
				xSeries[i] = new double[dimension];
			int index = 0;
			record(index, t0, x0);
		}

		public void record(double t, double[] x) {
			int index = (int)FastMath.ceil((numOfTimepoints - 1) * t / (t1 - t0));
			if (index >= numOfTimepoints)
				index = numOfTimepoints - 1;
			record(index, t, x);
		}

		private void record(int index, double t, double[] x) {
//			if (index < nextIndex)
//				nextIndex = index;
//			while (nextIndex < index) {
//				int prevIndex = nextIndex - 1;
//				for (int s=0; s < x.length; s++)
//					xSeries[lastIndexUsed][s] = x[s] * getScaleFactor(s);
//				nextIndex++;
//			}
			tSeries[index] = t;
			for (int i=0; i < dimension; i++)
				xSeries[index][i] = x[i];
			mask[index] = true;
//			nextIndex++;
		}

	};
//	private SimulationInformation simInfo;

	private class RootFunction implements UnivariateFunction {

    	private double offset = 0.0;
		private TemporaryTrajectoryInterpolator tti;
		private int index;

    	public RootFunction(TemporaryTrajectoryInterpolator tti, int index) {
    		this.tti = tti;
    		this.index = index;
    	}

    	public void setOffset(double offset) {
    		this.offset = offset;
    	}

    	public double getOffset() {
    		return offset;
    	}
   
		@Override
		public double value(double t) {
			double state = tti.getInterpolatedState(t, index);
			double v = state - offset;
			return v;
		}

		public void addOffset(double dOffset) {
			offset += dOffset;
		}

	}

	public double simulate(PDMPModel _model, final double t0, final double[] x0, double t1, double[] x1) {
		final int STOCHASTIC_CHECK_INTERVAL = 100;
		final int STOCHASTIC_COUNTDOWN_INTERVAL = 100;
		final double MIN_DETERMINISTIC_TO_COUPLED_RATIO = 2.0;
		final double STOP_TIME_ALPHA = 0.80;
//		rdg = new RandomDataGenerator();
//		rdg.reSeed(105L);

		AdaptiveMSHRNModel model = (AdaptiveMSHRNModel)_model;
		AdaptiveMSHRN net = (AdaptiveMSHRN)model.getNetwork();

		checkArgument(model.getNumberOfSpecies() == x0.length, "Expected model.getNumberOfSpecies() == x0.length but found %s != %s",
				model.getNumberOfSpecies(), x0.length);

		double t = t0;
		double[] x = new double[x0.length + 1];
		for (int i=0; i < x0.length; i++)
			x[i] = x0[i];
		double[] xDot = new double[x.length];
    	model.initialize(t, x);

        double[] xx = new double[x.length];

    	// Simulation information
		StateObserverAdapter stateObserver;
		boolean recordSimulationInformation = simulationInformationTrajectoryRecorders.size() > 0 || printMessages;
		stateObserver = new StateObserverAdapter(model, trajectoryRecorders, optionalTrajectoryRecorders);
//		this.simInfo = simInfo;
		stateObserver.initialize(t, x, t1);
		List<PDMPEventObserver> optionalEventObservers = model.getOptionalEventObservers();
		List<PDMPEventObserver> eventObservers = new ArrayList<PDMPEventObserver>(1 + optionalEventObservers.size());
//		eventObservers.add(model.getJumpEventObserver());
		eventObservers.addAll(optionalEventObservers);
		OdeAdapter ode = new OdeAdapter(model.getVectorField());
		EventObserverAdapter eventObserver = new EventObserverAdapter(eventObservers);
		eventObserver.initialize(t, x, t1);
		StochasticReactionNetworkModel rnm = model.getTransitionMeasure();
    	FirstOrderDifferentialEquations vectorField = model.getVectorField();

		double[] propVec = new double[rnm.getNumberOfReactions()];
    	double[] reactionCounterArray = new double[model.getNumberOfReactions()];
    	double msgDt = (t1 - t0) / 20.0;
    	double nextMsgT = t0 + msgDt;
    	int j = 0;
//    	CVodeSolver solver = (CVodeSolver)this.solver;
    	EulerSolver solver = (EulerSolver)this.solver;
    	TemporaryTrajectoryInterpolator tti = new TemporaryTrajectoryInterpolator(1001, x.length);
    	InnerTrajectoryRecorder itr = new InnerTrajectoryRecorder(1001);
    	double[] pm = new double[x.length];
    	model.computePrimaryState(t, x, pm);
    	itr.init(t0, pm, t1);
        double[] q = new double[x0.length];
        double[] p = new double[x0.length];

        double tLostTotal = 0.0;

		RootFunction rootFunction = new RootFunction(tti, x.length - 1);
//		UnivariateSolver rootSolver = new BrentSolver();
    	UnivariateSolver rootSolver = new BracketingNthOrderBrentSolver();

    	int stochasticCountdown = 0;

    	long ev = 0;
		synchronized (solver) {
			solver.initialize(ode, eventObserver, stateObserver, eventObserver);
			final long startTime = System.currentTimeMillis();
			while (true) {
				if (showProgress)
					while (t > nextMsgT) {
						System.out.println("Progress: " + (100 * (nextMsgT - t0)/(t1 - t0)) + "%");
						nextMsgT += msgDt;
					}
				boolean hasDeterministicPart = stochasticCountdown == 0 && model.hasVectorField();
				// TODO: Check whether this gives any performance gain
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
				double propSum = 0.0;
				boolean propensitiesComputed = false;
				double nextStopTime = t;

				if (hasDeterministicPart) {


					double deterministicPropensitiesSum = model.computeDeterministicPropensitiesSum(t, x);
//					double uncoupledPropensitiesSum = model.computeUncoupledPropensitiesSum(t, x);
					double coupledPropensitiesSum = model.computeCoupledPropensitiesSum(t, x);
//					double[] propensities = model.computeAllPropensities(t, x);

					if (deterministicPropensitiesSum / coupledPropensitiesSum < MIN_DETERMINISTIC_TO_COUPLED_RATIO) {
						hasDeterministicPart = false;
						stochasticCountdown = STOCHASTIC_COUNTDOWN_INTERVAL;
					} else {

			        	double integrationTimeStep;
			        	final double alpha = STOP_TIME_ALPHA;
			        	if (coupledPropensitiesSum > 0.0)
				        	integrationTimeStep = - FastMath.log(1 - alpha) / coupledPropensitiesSum;
		//	        		integrationTimeStep = 1 / coupledPropensitiesSum;
			        	else
			        		integrationTimeStep = 1 * (t1 - t);

			        	nextStopTime = FastMath.min(t + integrationTimeStep, t1);
					}
				}

		        if (hasDeterministicPart) {

					// Evolve ODE until next stochastic reaction fires
//		        	double nextUnitJumpTime = -FastMath.log(rdg.nextUniform(0.0,  1.0));
		        	// x[x.length - 1] will be < 0
//			        x[x.length - 1] = -nextUnitJumpTime;
		        	double tPrevious = t;
		        	x[x.length - 1] = 0.0;
			        model.checkAndHandleOptionalEvent(t, x);
		        	solver.prepareStep(t, x, t1);
		        	model.computePrimaryState(t, x, pm);
		        	pm[x.length - 1] = x[x.length - 1];
		        	tti.reset(t, pm);
		        	while (t < nextStopTime) {
//				        double[] pm1 = model.computePrimaryState(t, x).clone();
//			        	t = solver.integrate();
		        		ev = ode.getEvaluations();
//		        		double tStep = t + 1;
		        		double tStep = 10;
				        double tNext = solver.integrateStep(tStep);
//				        double[] pm2 = model.computePrimaryState(t, x).clone();
				        long newEv = ode.getEvaluations() - ev;

//				        long samples = FastMath.min(newEv, 10);
				        long samples = newEv;
				        double dt = (tNext - t) / (double)samples;
//				        double tt = t;
//				        for (int i=1; i < samples; i++) {
//				        	tt += dt;
//				        	solver.computeInterpolatedSolution(tt, xx);
//				        	model.computePrimaryState(tt, xx, pm);
//				        	pm[xx.length - 1] = xx[xx.length - 1];
//							tti.record(tt, pm);
//				        }

			        	model.computePrimaryState(t, x, pm);
			        	pm[x.length - 1] = x[x.length - 1];
				        tti.record(tNext, pm);
				    	itr.record(tNext, pm);
			        	t = tNext;
//				        solver.computeInterpolatedSolution(tHalf, xHalf);
//				        stateObserver.report(tHalf, xHalf);
//				        stateObserver.report(tNext, x);
//			        	if (showProgress)
//							while (t > nextMsgT) {
//								System.out.println("Progress: " + (100 * (nextMsgT - t0)/(t1 - t0)) + "%");
//								nextMsgT += msgDt;
//							}
			        	model.checkOptionalEvent(t, x);
			        	if (model.hasOptionalEventOccured()) {
				        	break;
//				        	solver.prepareOneStep(t, x);
			        	}
			        }

		        	// TODO
//			        System.out.println(String.format("Evolved from t=%f to t=%f", tPrevious, t));
			        // Find all reactions that fired between tPrevious and t
		        	int k = 0;
			        double tReaction = tPrevious;
		        	if (tReaction >= 2.6617299629380712)
		        		k = k;
			        double tt = t;
			        Arrays.fill(q, 0.0);
			        int coupledReactionCounter = 0;
			        int unCoupledReactionCounter = 0;
			        double offset = 0.0;
					while (tReaction < t) {
						double logU = -FastMath.log(rdg.nextUniform(0.0,  1.0));;
			        	double nextUnitJumpTime = logU;
			        	offset += nextUnitJumpTime;
			        	rootFunction.setOffset(offset);
//				        	final UnivariateInterpolator interpolator = new LinearInterpolator();
//				        	final UnivariateFunction propensityFunction = interpolator.interpolate(tSeries, xSeries[xSeries.length - 1]);
//				        	System.out.println(String.format("  rootFunction(tReaction)=%f, propensityFunction(tReaction)=%f, nextUnitJumpTime=%f", rootFunction.value(tReaction), propensityFunction.value(tReaction), nextUnitJumpTime));
//				        	System.out.println(String.format("  rootFunction(t)=%f, propensityFunction(t)=%f, nextUnitJumpTime=%f", rootFunction.value(t), propensityFunction.value(t), nextUnitJumpTime));
//				        	UnivariateFunction propensityFunction = new UnivariateFunction() {
	//
//				        		private int findIndex(double t) {
//				        			int index = Collections.binarySearch(recordList, new Record(t, new double[0]), new Comparator<Record>() {
	//
//										@Override
//										public int compare(Record o1, Record o2) {
//											return Double.compare(o1.t, o2.t);
//										}
	//
//				        			});
//				        			if (index >= 0)
//				        				return index;
//				        			else
//				        				return -index - 1;
//				        		}
	//
//								@Override
//								public double value(double t) {
//									int index = findIndex(t);
//									if (index + 1 < recordList.size()) {
//										double[] x1 = recordList.get(index).x;
//										double[] x2 = recordList.get(index + 1).x;
//										double t1 = recordList.get(index).t;
//										double t2 = recordList.get(index + 1).t;
//										double w1 = (t - t1) / (t2 - t1);
//										double w2 = 1.0 - t1;
//										double v = w1 * x1[x1.length - 1] + w2 * x2[x2.length - 2] - nextUnitJumpTime;
//										index++;
//										return v - nextUnitJumpTime;
//									} else if (index < recordList.size()) {
//										double[] x1 = recordList.get(index).x;
//										double v = x1[x1.length - 1];
//										return v - nextUnitJumpTime;
//									} else
//										return 0.0;
//								}
//			
//							};
			        	try {
							if (rootFunction.value(t) < 0) {
//									System.out.println(String.format("  No stochastic reaction took place: rootFunction(t)=%f", rootFunction.value(t)));
								break;
							}
			        	} catch (Exception e) {
			        		rootFunction.value(t);
			        	}

			        	try {
			        		tReaction = rootSolver.solve(1000, rootFunction, tReaction, t);
//			        		System.out.println(String.format("  found root at tReaction=%f, rootFunction(tReaction)=%f, propensityFunction(tReaction)=%f", tReaction, rootFunction.value(tReaction), propensityFunction.value(tReaction)));
			        	} catch (Exception e) {
//				        	System.out.println(String.format("  rootFunction(tReaction)=%f, propensityFunction(tReaction)=%f, nextUnitJumpTime=%f", rootFunction.value(tReaction), propensityFunction.value(tReaction), nextUnitJumpTime));
//				        	System.out.println(String.format("  rootFunction(t)=%f, propensityFunction(t)=%f, nextUnitJumpTime=%f", rootFunction.value(t), propensityFunction.value(t), nextUnitJumpTime));
		        			tReaction = rootSolver.solve(1000, rootFunction, tReaction, t);
			        	}
			        	double[] tmp = tti.getInterpolatedState(tReaction);
//			        	computeInterpolatedScaledState(tReaction, x, net);
			        	for (int i=0; i < q.length; i++)
			        		p[i] = (tmp[i] + q[i]) * net.getInverseSpeciesScaleFactor(i);
				        propSum = model.computePropensitiesAndSum(tReaction, p, propVec);
//				        rnm.computePropensities(tReaction, p, propVec);
				        int reaction = DiscreteProbabilityDistribution.sample(rdg, propVec, propSum);
//							System.out.println(String.format("  Stochastic reaction r=%d took place at t=%f", reaction, tReaction));
				        if (model.coupledStochasticReactions.contains(reaction)) {
//				        if (coupledStochasticReactionMask[reaction]) {
//					        	System.out.println("  Stochastic reaction was coupled to vector field");
				        	t = tReaction;
//							System.out.println(String.format("  Rewinding to t=%f", t));
				        	coupledReactionCounter++;
				        } else {
				        	unCoupledReactionCounter++;
				        }
				        rnm.changeState(reaction, t, q);
			        	reactionCounterArray[reaction]++;
			        	if (k > STOCHASTIC_CHECK_INTERVAL) {
			        		model.computePrimaryState(tReaction, p, pm);
				        	itr.record(tReaction, pm);
			        		model.checkAndHandleOptionalEvent(tReaction, p);
			        		model.manualUpdate(propVec);
			        		k = 0;
			        	}
			        	k++;

			        	// Update the integrated sum of propensities
	        			int index = tti.binarySearch(tReaction);
						if (index < 0)
							index = (-index - 1) - 1;
						index += 1;
						double timepoint = tReaction;
						double[] state = tti.getInterpolatedState(tReaction);
						double y = state[x.length - 1];
						for (int i=index; i < tti.getNumberOfRecordedTimepoints(); i++) {
							double dt = tti.getRecordedTimepoint(index) - timepoint;
				        	for (int l=0; l < p.length; l++)
				        		p[l] = (state[l] + q[l]) * net.getInverseSpeciesScaleFactor(l);
							double f = model.computePropensitySum(tti.getRecordedTimepoint(index), p);
							double dy = f * dt;
							y += dy;
							tti.getRecordedState(index)[x.length - 1] = y;
							timepoint = tti.getRecordedTimepoint(index);
							state = tti.getRecordedState(index);
						}

						// Update offset for finding reaction firing time
		        		offset += rootFunction.value(tReaction);
					}

//					System.out.println(String.format("  total=%d, coupled=%d, uncoupled=%d", coupledReactionCounter + unCoupledReactionCounter, coupledReactionCounter, unCoupledReactionCounter));
					double tLost = tt - tReaction;
					tLostTotal += tLost;
					double dt = tt - tPrevious;
//					System.out.println(String.format("Evolved from %f to %f but coupled reaction happened at %f [lost %f of %f, %f%%]", tPrevious, tt, tReaction, tLost, dt, 100.0 * tLost / dt));

//					computeInterpolatedScaledState(t, x, net);
//		        	for (int i=0; i < q.length; i++)
//		        		x[i] = tmp[i] + q[i];

					if (coupledReactionCounter > 0) {
						double[] tmp = tti.getInterpolatedState(t);
			        	for (int i=0; i < q.length; i++) {
			        		tmp[i] += q[i];
			        		x[i] = tmp[i] * net.getInverseSpeciesScaleFactor(i);
			        	}
						itr.record(t, tmp);
					} else if (unCoupledReactionCounter > 0) {
			        	for (int i=0; i < q.length; i++) {
			        		x[i] += q[i] * net.getInverseSpeciesScaleFactor(i);
			        	}
			        	model.computePrimaryState(t, x, pm);
						itr.record(t, pm);
					}

		        	if (model.hasOptionalEventOccured())
			        	model.handleOptionalEvent(t, x);

				} else {

			        propSum = model.computeAllPropensitiesAndSum(t, x, propVec);
			        // Find next reaction time point
			        if (propSum < 0)
			        	throw new RuntimeException("TODO");
			        if (propSum == 0) {
			        	t = t1;
			        	break;
			        }

			        double tau = rdg.nextExponential(1 / propSum);
			        t = t + tau;

			        int reaction = DiscreteProbabilityDistribution.sample(rdg, propVec, propSum);

			        if (reaction >= 0) {
			    		rnm.changeState(reaction, t, x);
			        	reactionCounterArray[reaction]++;
//			        	stateObserver.report(t, x);
			        	// TODO: Make this value configurable (coupled to N?)
			        	if (j > STOCHASTIC_CHECK_INTERVAL) {
			        		model.computePrimaryState(t, x, pm);
				        	itr.record(t, pm);
			        		model.checkAndHandleOptionalEvent(t, x);
			        		model.manualUpdate(propVec);
			        		j = 0;
			        	}
			        	j++;
			        } else
			        	throw new RuntimeException("TODO");

					if (stochasticCountdown > 0)
						stochasticCountdown--;
				}
	
		        // Stop if we reached the end-timepoint
		        if (t >= t1)
		        	break;

			}
			if (printMessages) {
				System.out.println("Integration time lost in total: " + tLostTotal);
//				System.out.println("Integrator invocations: " + simInfo.getIntegrationSteps());
				System.out.println("Observed " + eventObserver.getEventCount() + " events");
				final long endTime = System.currentTimeMillis();
				System.out.println("Execution time: " + (endTime - startTime));
		    	long evaluationCounter = ode.getEvaluations();
//				System.out.println("Total of " + evaluationCounter + " evaluations and " + simInfo.getReactionCount() + " reactions performed");
				Utilities.printArray("Total reaction counts", reactionCounterArray);
				for (int r=0; r < reactionCounterArray.length; r++)
//					reactionCounterArray[r] /= simInfo.getReactionCount();
				Utilities.printArray("Relative reaction counts", reactionCounterArray);
			}
			for (int i=0; i < x1.length; i++)
				x1[i] = x[i];
//	    	for (TrajectoryRecorder<T> handler : trajectoryRecorders)
//	    		handler.setFinalState(t1, x1);
//        	stateObserver.report(t1, x1);
    		model.computePrimaryState(t, x, pm);
			itr.record(t, pm);
	    	solver.dispose();

	    	double[][] qSeries = new double[itr.xSeries[0].length][];
	    	for (int i=0; i < qSeries.length; i++) {
	    		qSeries[i] = new double[itr.tSeries.length];
	    		for (int ii=0; ii < itr.tSeries.length; ii++)
	    			qSeries[i][ii] = itr.xSeries[ii][i];
	    	}

	    	double it = itr.tSeries[0];
	    	double[] ix = itr.xSeries[0];
	    	for (TrajectoryRecorder tr : trajectoryRecorders) {
	    		tr.beginRecording(t0, x0, t1);
	    		for (int i=0; i < itr.tSeries.length; i++) {
	    			boolean b = itr.mask[i];
	    			if (b) {
	    				if (itr.xSeries[i][2] < ix[2])
	    					t = t;
		    			it = itr.tSeries[i];
		    			ix = itr.xSeries[i];
	    			}
	    			tr.record(it, ix);
	    		}
	    		tr.endRecording(x1);
	    	}
		}
		return t;
	}

	@Override
	public void addTrajectoryRecorder(TrajectoryRecorder tr) {
		trajectoryRecorders.add(tr);
	}

	@Override
	public void removeTrajectoryRecorder(TrajectoryRecorder tr) {
		trajectoryRecorders.remove(tr);
	}

	@Override
	public void clearTrajectoryRecorders() {
		trajectoryRecorders.clear();
	}

	public void addOptionalTrajectoryRecorder(TrajectoryRecorder tr) {
		optionalTrajectoryRecorders.add(tr);
	}

	public void removeOptionalTrajectoryRecorder(TrajectoryRecorder tr) {
		optionalTrajectoryRecorders.remove(tr);
	}

	public void clearOptionalTrajectoryRecorders() {
		optionalTrajectoryRecorders.clear();
	}

	public void addSimulationInformationTrajectoryRecorder(TrajectoryRecorder tr) {
		simulationInformationTrajectoryRecorders.add(tr);
	}

	public void removeSimulationInformationTrajectoryRecorder(TrajectoryRecorder tr) {
		simulationInformationTrajectoryRecorders.remove(tr);
	}

	public void clearSimulationInformationTrajectoryRecorders() {
		simulationInformationTrajectoryRecorders.clear();
	}

}
