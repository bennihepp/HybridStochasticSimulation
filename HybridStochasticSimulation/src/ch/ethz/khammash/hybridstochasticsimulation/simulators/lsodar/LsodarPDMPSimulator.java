package ch.ethz.khammash.hybridstochasticsimulation.simulators.lsodar;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.util.FastMath;

import ch.ethz.khammash.hybridstochasticsimulation.Utilities;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPEventObserver;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.ContinuousTrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;
import ch.ethz.khammash.nativeode.Solver;
import ch.ethz.khammash.nativeode.lsodar.LsodarSolver;

public class LsodarPDMPSimulator<T extends PDMPModel> implements Simulator<T, ContinuousTrajectoryRecorder<T>> {

	public static final double DEFAULT_EVENT_HANDLER_MAX_CHECK_INTERVAL = Double.POSITIVE_INFINITY;
	public static final double DEFAULT_EVENT_HANDLER_CONVERGENCE = 1e-12;
	public static final double DEFAULT_EVENT_HANDLER_CONVERGENCE_FACTOR = Double.NaN;
	public static final int DEFAULT_EVENT_HANDLER_MAX_ITERATION_COUNT = 1000;

	private RandomDataGenerator rdg;
	private List<FiniteTrajectoryRecorder<T>> trajectoryRecorders;
	private Solver solver;

	public LsodarPDMPSimulator() {
		this(null);
	}

	public LsodarPDMPSimulator(RandomDataGenerator rdg) {
		if (rdg == null)
			rdg = new RandomDataGenerator();
		solver = LsodarSolver.getInstance();
		this.rdg = rdg;
		trajectoryRecorders = new LinkedList<FiniteTrajectoryRecorder<T>>();
	}

	public double simulate(T model, final double t0, final double[] x0, double t1, double[] x1) {
		rdg = new RandomDataGenerator();
		rdg.reSeed(10L);
		boolean showProgress = false;
		double t = t0;
		synchronized (solver) {
			for (FiniteTrajectoryRecorder<T> tr : trajectoryRecorders)
				tr.setModel(model);
			LsodarStateObserverAdapter<T> stateObserver = new LsodarStateObserverAdapter<T>(trajectoryRecorders);
			stateObserver.initialize(t0, x0, t1);
			List<PDMPEventObserver> optionalEventObservers = model.getOptionalEventObservers();
			List<PDMPEventObserver> eventObservers = new ArrayList<PDMPEventObserver>(1 + optionalEventObservers.size());
			eventObservers.add(model.getJumpEventObserver());
			eventObservers.addAll(optionalEventObservers);
			LsodarOdeAdapter<T> ode = new LsodarOdeAdapter<T>(model);
			LsodarEventObserverAdapter eventObserver = new LsodarEventObserverAdapter(eventObservers);
			eventObserver.initialize(t0, x0, t1);
			solver.initialize(ode, eventObserver, stateObserver, eventObserver);
	    	StochasticReactionNetworkModel rnm = model.getTransitionMeasure();
	    	FirstOrderDifferentialEquations vectorField = model.getVectorField();
	    	LsodarTimepointIterator<T> timepointIterator;
			if (trajectoryRecorders.size() > 0)
	    		timepointIterator = new LsodarTimepointIterator<T>(trajectoryRecorders.iterator().next());
	    	else
	    		timepointIterator = new LsodarTimepointIterator<T>(t0, t1);
			double[] x = new double[x0.length + 2];
			for (int i=0; i < x0.length; i++)
				x[i] = x0[i];
			double[] xDot = x.clone();
	    	model.initialize(t, x);
			double[] propVec = new double[rnm.getNumberOfReactions()];
	    	long reactionCounter = 0;
	    	double[] reactionCounterArray = new double[model.getNumberOfReactions()];
	    	double msgDt = (t1 - t0) / 20.0;
	    	double nextMsgT = t0 + msgDt;
//	    	int j = 0;
			final long startTime = System.currentTimeMillis();
			while (t < t1) {
				if (showProgress)
					while (t > nextMsgT) {
						System.out.println("Progress: " + (100 * (nextMsgT - t0)/(t1 - t0)) + "%");
						nextMsgT += msgDt;
					}
				boolean hasDeterministicPart = model.hasVectorField();
				if (hasDeterministicPart && model.isTimeIndependent()) {
					vectorField.computeDerivatives(t, x, xDot);
					boolean allZero = true;
					for (int s=0; s < x0.length; s++)
						if (xDot[s] != 0.0) {
							allZero = false;
							break;
						}
						if (allZero)
							hasDeterministicPart = false;
				}
				double propSum = 0.0;
				boolean propensitiesComputed = false;
				if (hasDeterministicPart) {
					// Evolve ODE until next stochastic reaction fires
			        x[x.length - 2] = 0.0;
			        x[x.length - 1] = -FastMath.log(rdg.nextUniform(0.0,  1.0));
//			        do {
//			        	model.handleOptionalEvent(t, x);
			        	timepointIterator.setCurrentTime(t);
			        	t = solver.integrate(timepointIterator, x);
			        	if (showProgress)
							while (t > nextMsgT) {
								System.out.println("Progress: " + (100 * (nextMsgT - t0)/(t1 - t0)) + "%");
								nextMsgT += msgDt;
							}
//			        } while (model.hasOptionalEventOccured());
				} else {
			        rnm.computePropensities(t, x, propVec);
			        for (int i=0; i < propVec.length; i++)
			        	propSum += propVec[i];
			        propensitiesComputed = true;
			        // Find next reaction time point
			        if (propSum <= 0)
			        	break;
			        // -FastMath.log(rng.nextUniform(0.0,  1.0))
			        double tau = rdg.nextExponential(1 / propSum);
			        t = t + tau;
				}
	
		        // Stop if we reached the end-timepoint
		        if (t >= t1)
		        	break;
	
		        if (!propensitiesComputed) {
			        // Determine which reaction fired and update state
			        rnm.computePropensities(t, x, propVec);
			        for (int i=0; i < propVec.length; i++)
			        	propSum += propVec[i];
		        }
		        double u = rdg.nextUniform(0.0, 1.0);
		        double w = 0.0;
		        int reaction = -1;
		        for (int l=0; l < propVec.length; l++) {
		        	w = w + propVec[l] / propSum;
		        	if (u < w) {
		        		reaction = l;
		        		rnm.updateState(reaction, t, x);
		        		break;
		        	}
		        }
		        if (reaction >= 0) {
		        	reactionCounter++;
		        	reactionCounterArray[reaction]++;
	//	        	for (TrajectoryRecorder<T> handler : trajectoryRecorders)
	//	        		handler.reportState(t, x);
		        	stateObserver.report(t, x);
//		        	if (j > 1000)
//		        		model.checkAndHandleOptionalEvent(t, x);
//		        	j++;
		        }
			}
			System.out.println("Observed " + eventObserver.getEventCount() + " events");
			final long endTime = System.currentTimeMillis();
			System.out.println("Execution time: " + (endTime - startTime));
	    	long evaluationCounter = ode.getEvaluations();
			System.out.println("Total of " + evaluationCounter + " evaluations and " + reactionCounter + " reactions performed");
			Utilities.printArray("Total reaction counts", reactionCounterArray);
			for (int r=0; r < reactionCounterArray.length; r++)
				reactionCounterArray[r] /= reactionCounter;
			Utilities.printArray("Relative reaction counts", reactionCounterArray);
			for (int i=0; i < x1.length; i++)
				x1[i] = x[i];
	    	for (TrajectoryRecorder<T> handler : trajectoryRecorders)
	    		handler.setFinalState(t1, x1);
	    	solver.dispose();
		}
		return t;
	}

	@Override
	public void addTrajectoryRecorder(ContinuousTrajectoryRecorder<T> tr) {
		if (tr instanceof FiniteTrajectoryRecorder<?>) {
			@SuppressWarnings("unchecked")
			FiniteTrajectoryRecorder<T> ftr = (FiniteTrajectoryRecorder<T>)tr;
			trajectoryRecorders.add(ftr);
		}
	}

	@Override
	public void removeTrajectoryRecorder(ContinuousTrajectoryRecorder<T> tr) {
		trajectoryRecorders.remove(tr);
	}

	@Override
	public void clearTrajectoryRecorders() {
		trajectoryRecorders.clear();
	}

}
