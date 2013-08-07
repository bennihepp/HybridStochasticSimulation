package ch.ethz.khammash.hybridstochasticsimulation.simulators;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.util.FastMath;

import ch.ethz.khammash.hybridstochasticsimulation.Utilities;
import ch.ethz.khammash.hybridstochasticsimulation.math.RandomDataUtilities;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.ode.EventObserverAdapter;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.ode.ExtendedStateObserverAdapter;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.ode.OdeAdapter;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.ode.StateObserverAdapter;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;
import ch.ethz.khammash.ode.FiniteTimepointProvider;
import ch.ethz.khammash.ode.Solver;
import ch.ethz.khammash.ode.lsodar.LsodarDirectSolver;

public class PDMPSimulator extends AbstractSimulator<PDMPModel> {

	public static final double DEFAULT_EVENT_HANDLER_MAX_CHECK_INTERVAL = Double.POSITIVE_INFINITY;
	public static final double DEFAULT_EVENT_HANDLER_CONVERGENCE = 1e-12;
	public static final double DEFAULT_EVENT_HANDLER_CONVERGENCE_FACTOR = Double.NaN;
	public static final int DEFAULT_EVENT_HANDLER_MAX_ITERATION_COUNT = 1000;

	private RandomDataGenerator rdg;
	private List<TrajectoryRecorder> trajectoryRecorders;
	private List<TrajectoryRecorder> optionalTrajectoryRecorders;
	private List<TrajectoryRecorder> simulationInformationTrajectoryRecorders;
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
		optionalTrajectoryRecorders = new LinkedList<TrajectoryRecorder>();
		simulationInformationTrajectoryRecorders = new LinkedList<TrajectoryRecorder>();
	}

	private class DefaultSimulationInformation implements SimulationInformation {
    	private long integrationCount = 0;
    	private long reactionCount = 0;
    	private boolean integrating = false;
		@Override
		public long getIntegrationCount() {
			return integrationCount;
		}
		@Override
		public long getReactionCount() {
			return reactionCount;
		}
		@Override
		public boolean isIntegrating() {
			return integrating;
		}
		@Override
		public double[] computeInformationState() {
			double[] state = new double[3];
			state[0] = integrating ? 1.0 : -1.0;
			state[1] = integrationCount;
			state[2] = reactionCount;
			return state;
		}
	}

	public double simulate(PDMPModel model, final double t0, final double[] x0, double t1, double[] x1) {
//		rdg = new RandomDataGenerator();
//		rdg.reSeed(105L);

		double t = t0;
		double[] x = new double[x0.length + 1];
		for (int i=0; i < x0.length; i++)
			x[i] = x0[i];
		double[] xDot = new double[x.length];
    	model.initialize(t, x);

    	// Simulation information
		StateObserverAdapter stateObserver;
		boolean recordSimulationInformation = simulationInformationTrajectoryRecorders.size() > 0 || printMessages;
    	DefaultSimulationInformation simInfo = null;
		if (recordSimulationInformation) {
	    	simInfo = new DefaultSimulationInformation();
//	    	simInfo.integrationCount = 0;
//	    	simInfo.reactionCount = 0;
//	    	simInfo.integrating = false;
			stateObserver = new ExtendedStateObserverAdapter(simInfo, model, trajectoryRecorders, optionalTrajectoryRecorders, simulationInformationTrajectoryRecorders);
		} else
			stateObserver = new StateObserverAdapter(model, trajectoryRecorders, optionalTrajectoryRecorders);
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
					timepointProvider = new FiniteTimepointProvider(ftr.gettSeries());
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
		synchronized (solver) {
			solver.initialize(ode, eventObserver, stateObserver, eventObserver);
			final long startTime = System.currentTimeMillis();
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
		        	if (recordSimulationInformation)
		        		simInfo.integrating = true;
					// Evolve ODE until next stochastic reaction fires
		        	double nextUnitJumpTime = -FastMath.log(rdg.nextUniform(0.0,  1.0));
		        	// x[x.length - 1] will be < 0
			        x[x.length - 1] = -nextUnitJumpTime;
			        model.checkAndHandleOptionalEvent(t, x);
		        	timepointProvider.setCurrentTimepoint(t);
		        	solver.prepare(timepointProvider, x);
		        	while (true) {
			        	t = solver.integrate();
			        	if (recordSimulationInformation)
			        		simInfo.integrationCount++;
			        	if (showProgress)
							while (t > nextMsgT) {
								System.out.println("Progress: " + (100 * (nextMsgT - t0)/(t1 - t0)) + "%");
								nextMsgT += msgDt;
							}
			        	if (model.hasOptionalEventOccured()) {
				        	model.handleOptionalEvent(t, x);
				        	timepointProvider.setCurrentTimepoint(t);
				        	solver.prepare(timepointProvider, x);
			        	} else
			        		break;
			        }
		        	if (recordSimulationInformation)
		        		simInfo.integrating = false;
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

		        // Determine which reaction fired and update state
		        if (!propensitiesComputed) {
			        rnm.computePropensities(t, x, propVec);
			        for (int i=0; i < propVec.length; i++)
			        	propSum += propVec[i];
		        }
		        int reaction = RandomDataUtilities.sampleFromProbabilityMassFunction(rdg, propVec);
	    		rnm.changeState(reaction, t, x);
//		        double u = rdg.nextUniform(0.0, 1.0);
//		        double w = 0.0;
//		        int reaction = -1;
//		        for (int l=0; l < propVec.length; l++) {
//		        	w = w + propVec[l] / propSum;
//		        	if (u < w) {
//		        		reaction = l;
//		        		rnm.changeState(reaction, t, x);
//		        		break;
//		        	}
//		        }
		        if (reaction >= 0) {
		        	if (recordSimulationInformation)
		        		simInfo.reactionCount++;
		        	reactionCounterArray[reaction]++;
		        	stateObserver.report(t, x);
		        	// TODO: Make this value configurable (coupled to N?)
		        	if (j > 100) {
		        		model.checkAndHandleOptionalEvent(t, x);
		        		j = 0;
		        	}
		        	j++;
		        }
			}
			if (printMessages) {
				System.out.println("Integrator invocations: " + simInfo.getIntegrationCount());
				System.out.println("Observed " + eventObserver.getEventCount() + " events");
				final long endTime = System.currentTimeMillis();
				System.out.println("Execution time: " + (endTime - startTime));
		    	long evaluationCounter = ode.getEvaluations();
				System.out.println("Total of " + evaluationCounter + " evaluations and " + simInfo.getReactionCount() + " reactions performed");
				Utilities.printArray("Total reaction counts", reactionCounterArray);
				for (int r=0; r < reactionCounterArray.length; r++)
					reactionCounterArray[r] /= simInfo.getReactionCount();
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
