package ch.ethz.khammash.hybridstochasticsimulation.simulators;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.util.FastMath;

import ch.ethz.khammash.hybridstochasticsimulation.Utilities;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.ode.EventObserverAdapter;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.ode.FiniteTimepointProvider;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.ode.OdeAdapter;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.ode.StateObserverAdapter;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;
import ch.ethz.khammash.ode.Solver;
import ch.ethz.khammash.ode.lsodar.LsodarDirectSolver;

public class PDMPSimulator implements Simulator<PDMPModel> {

	public static final double DEFAULT_EVENT_HANDLER_MAX_CHECK_INTERVAL = Double.POSITIVE_INFINITY;
	public static final double DEFAULT_EVENT_HANDLER_CONVERGENCE = 1e-12;
	public static final double DEFAULT_EVENT_HANDLER_CONVERGENCE_FACTOR = Double.NaN;
	public static final int DEFAULT_EVENT_HANDLER_MAX_ITERATION_COUNT = 1000;

	private RandomDataGenerator rdg;
	private List<TrajectoryRecorder> trajectoryRecorders;
	private Solver solver;

	public PDMPSimulator() {
		this(null, null);
	}

	public PDMPSimulator(RandomDataGenerator rdg) {
		this(null, rdg);
	}

	public PDMPSimulator(Solver solver) {
		this(solver, null);
	}

	public PDMPSimulator(Solver solver, RandomDataGenerator rdg) {
		if (rdg == null)
			rdg = new RandomDataGenerator();
		this.rdg = rdg;
		if (solver == null)
			solver = LsodarDirectSolver.getInstance();
		this.solver = solver;
		trajectoryRecorders = new LinkedList<TrajectoryRecorder>();
	}

	public int isIntegrating;
	public int integratorCounter;
	public int reactionCounter;

	public double simulate(PDMPModel model, final double t0, final double[] x0, double t1, double[] x1) {
//		rdg = new RandomDataGenerator();
//		rdg.reSeed(102L);
		boolean printMessages = false;
		boolean showProgress = false;

		double t = t0;
		double[] x = new double[x0.length + 2];
		for (int i=0; i < x0.length; i++)
			x[i] = x0[i];
		double[] xDot = new double[x.length];
    	model.initialize(t, x);

		StateObserverAdapter stateObserver = new StateObserverAdapter(model, trajectoryRecorders);
		stateObserver.initialize(t, x, t1);
		List<PDMPEventObserver> optionalEventObservers = model.getOptionalEventObservers();
		List<PDMPEventObserver> eventObservers = new ArrayList<PDMPEventObserver>(1 + optionalEventObservers.size());
		eventObservers.add(model.getJumpEventObserver());
		eventObservers.addAll(optionalEventObservers);
		OdeAdapter ode = new OdeAdapter(model);
		EventObserverAdapter eventObserver = new EventObserverAdapter(eventObservers);
		eventObserver.initialize(t, x, t1);
		StochasticReactionNetworkModel rnm = model.getTransitionMeasure();
    	FirstOrderDifferentialEquations vectorField = model.getVectorField();

    	FiniteTimepointProvider timepointProvider = null;
		if (trajectoryRecorders.size() > 0) {
			for (TrajectoryRecorder tr : trajectoryRecorders)
				if (tr instanceof FiniteTrajectoryRecorder) {
					FiniteTrajectoryRecorder ftr = (FiniteTrajectoryRecorder)tr; 
					timepointProvider = new FiniteTimepointProvider(ftr);
					break;
				}
		}
		if (timepointProvider == null)
			timepointProvider = new FiniteTimepointProvider(t0, t1);

		double[] propVec = new double[rnm.getNumberOfReactions()];
    	double[] reactionCounterArray = new double[model.getNumberOfReactions()];
    	double msgDt = (t1 - t0) / 20.0;
    	double nextMsgT = t0 + msgDt;
    	int j = 0;
    	integratorCounter = 0;
    	reactionCounter = 0;
    	isIntegrating = 0;
		synchronized (solver) {
			solver.initialize(ode, eventObserver, stateObserver, eventObserver);
			final long startTime = System.currentTimeMillis();
//			while (t < t1) {
			while (true) {
				if (showProgress)
					while (t > nextMsgT) {
						System.out.println("Progress: " + (100 * (nextMsgT - t0)/(t1 - t0)) + "%");
						nextMsgT += msgDt;
					}
				boolean hasDeterministicPart = model.hasVectorField();
				// TODO: Check whether this gives any performance gain
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
	        		isIntegrating = 1;
					// Evolve ODE until next stochastic reaction fires
			        x[x.length - 2] = 0.0;
			        x[x.length - 1] = -FastMath.log(rdg.nextUniform(0.0,  1.0));
			        model.checkAndHandleOptionalEvent(t, x);
		        	timepointProvider.setCurrentTime(t);
		        	solver.prepare(timepointProvider, x);
		        	while (true) {
			        	t = solver.integrate();
						integratorCounter++;
			        	if (showProgress)
							while (t > nextMsgT) {
								System.out.println("Progress: " + (100 * (nextMsgT - t0)/(t1 - t0)) + "%");
								nextMsgT += msgDt;
							}
			        	if (model.hasOptionalEventOccured()) {
				        	model.handleOptionalEvent(t, x);
				        	timepointProvider.setCurrentTime(t);
				        	solver.prepare(timepointProvider, x);
			        	} else
			        		break;
			        }
	        		isIntegrating = 0;
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
		        		rnm.changeState(reaction, t, x);
		        		break;
		        	}
		        }
		        if (reaction >= 0) {
		        	reactionCounter++;
		        	reactionCounterArray[reaction]++;
		        	stateObserver.report(t, x);
		        	// TODO: Should this also be coupled to N?
		        	if (j > 100) {
		        		model.checkAndHandleOptionalEvent(t, x);
		        		j = 0;
		        	}
		        	j++;
		        }
			}
			if (printMessages) {
				System.out.println("Integrator invocations: " + integratorCounter);
				System.out.println("Observed " + eventObserver.getEventCount() + " events");
				final long endTime = System.currentTimeMillis();
				System.out.println("Execution time: " + (endTime - startTime));
		    	long evaluationCounter = ode.getEvaluations();
				System.out.println("Total of " + evaluationCounter + " evaluations and " + reactionCounter + " reactions performed");
				Utilities.printArray("Total reaction counts", reactionCounterArray);
				for (int r=0; r < reactionCounterArray.length; r++)
					reactionCounterArray[r] /= reactionCounter;
				Utilities.printArray("Relative reaction counts", reactionCounterArray);
			}
			for (int i=0; i < x1.length; i++)
				x1[i] = x[i];
//	    	for (TrajectoryRecorder<T> handler : trajectoryRecorders)
//	    		handler.setFinalState(t1, x1);
        	stateObserver.report(t, x1);
	    	solver.dispose();
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

}
