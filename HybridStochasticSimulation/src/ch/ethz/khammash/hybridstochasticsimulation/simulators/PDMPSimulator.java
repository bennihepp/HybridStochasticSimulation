package ch.ethz.khammash.hybridstochasticsimulation.simulators;

import java.util.LinkedHashSet;
import java.util.Set;

import org.apache.commons.math3.ode.AbstractIntegrator;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.events.EventHandler;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.util.FastMath;

import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.ContinuousTrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;


public class PDMPSimulator<T extends PDMPModel> implements Simulator<T, ContinuousTrajectoryRecorder<T>> {

	public static final double DEFAULT_EVENT_HANDLER_MAX_CHECK_INTERVAL = Double.POSITIVE_INFINITY;
	public static final double DEFAULT_EVENT_HANDLER_CONVERGENCE = 1e-12;
	public static final double DEFAULT_EVENT_HANDLER_CONVERGENCE_FACTOR = Double.NaN;
	public static final int DEFAULT_EVENT_HANDLER_MAX_ITERATION_COUNT = 1000;

	protected RandomDataGenerator rdg;
	protected Set<ContinuousTrajectoryRecorder<T>> trajectoryRecorders;
	protected AbstractIntegrator integrator;
	private double ehMaxCheckInterval;
	private double ehConvergence;
	private double ehConvergenceFactor;
	private int ehMaxIterationCount;

	public PDMPSimulator() {
		this(null, null);
	}

	public PDMPSimulator(RandomDataGenerator rdg) {
		this(null, rdg);
	}

	public PDMPSimulator(AbstractIntegrator integrator) {
		this(integrator, null);
	}

	public PDMPSimulator(AbstractIntegrator integrator, RandomDataGenerator rdg) {
		ehMaxCheckInterval = DEFAULT_EVENT_HANDLER_MAX_CHECK_INTERVAL;
		ehConvergence = DEFAULT_EVENT_HANDLER_CONVERGENCE;
		ehConvergenceFactor = DEFAULT_EVENT_HANDLER_CONVERGENCE_FACTOR;
		ehMaxIterationCount = DEFAULT_EVENT_HANDLER_MAX_ITERATION_COUNT;
		if (integrator == null)
			integrator = new DormandPrince853Integrator(1.0e-8, 100.0, 1.0e-10, 1.0e-10);
		this.integrator = integrator;
		if (rdg == null)
			rdg = new RandomDataGenerator();
		this.rdg = rdg;
		trajectoryRecorders = new LinkedHashSet<>();
	}

	public double getEventHandlerMaxCheckInterval() {
		return ehMaxCheckInterval;
	}

	public void setEventHandlerMaxCheckInterval(double maxCheckInterval) {
		this.ehMaxCheckInterval = maxCheckInterval;
	}

	public double getEventHandlerConvergence() {
		return ehConvergence;
	}

	public void setEventHandlerConvergence(double convergence) {
		this.ehConvergence = convergence;
		if (!Double.isNaN(convergence))
			this.ehConvergenceFactor = Double.NaN;
	}

	public double getEventHandlerConvergenceFactor() {
		return ehConvergenceFactor;
	}

	public void setEventHandlerConvergenceFactor(double convergenceFactor) {
		this.ehConvergenceFactor = convergenceFactor;
		if (!Double.isNaN(convergenceFactor))
			this.ehConvergence = Double.NaN;
	}

	public int getEventHandlerMaxIterationCount() {
		return ehMaxIterationCount;
	}

	public void setEventHandlerMaxIterationCount(int maxIterationCount) {
		this.ehMaxIterationCount = maxIterationCount;
	}

	public double simulate(T model, final double t0, final double[] x0, double t1, double[] x1) {
		integrator.clearEventHandlers();
		integrator.clearStepHandlers();
//    	ExpandableStatefulODE ode = new ExpandableStatefulODE(model.getFirstOrderDifferentialEquations());
		FirstOrderDifferentialEquations ode = model.getVectorField();
    	StochasticReactionNetworkModel rnm = model.getTransitionMeasure();
    	EventHandler pdmpEventHandler = model.getJumpEventObserver();
		double[] x = new double[x0.length + 2];
		for (int i=0; i < x0.length; i++)
			x[i] = x0[i];
		double[] xDot = x.clone();
		double t = t0;
    	model.initialize(t, x);
		double[] propVec = new double[rnm.getNumberOfReactions()];
		for (ContinuousTrajectoryRecorder<T> tr : trajectoryRecorders) {
			tr.setModel(model);
			tr.setInitialState(t, x);
			tr.init(t0, x, t1);
			integrator.addStepHandler(tr);
		}
		double conv = Double.isNaN(ehConvergenceFactor) ? ehConvergence : ehConvergenceFactor * (t1 - t0);
		if (Double.isNaN(conv))
			throw new IllegalArgumentException("Either convergence or convergence factor must be a positive real number");
		integrator.addEventHandler(pdmpEventHandler, ehMaxCheckInterval, conv, ehMaxIterationCount);
    	for (EventHandler eh : model.getOptionalEventObservers())
    		integrator.addEventHandler(eh, ehMaxCheckInterval, conv, ehMaxIterationCount);
//    	long evaluationCounter = 0;
//    	long reactionCounter = 0;
//    	double[] reactionCounterArray = new double[model.getPropensityDimension()];
		while (t < t1) {
			boolean hasDeterministicPart = model.hasVectorField();
			if (hasDeterministicPart && model.isTimeIndependent()) {
				ode.computeDerivatives(t, x, xDot);
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
		        do {
		        	model.handleOptionalEvent(t, x);
		        	t = integrator.integrate(ode, t, x, t1, x);
//		        	ode.setTime(t);
//		        	ode.setPrimaryState(x);
//		        	integrator.integrate(ode, t1);
//		        	System.arraycopy(ode.getPrimaryState(), 0, x, 0, x.length);
//		        	t = ode.getTime();
//		        	evaluationCounter += integrator.getEvaluations();
		        } while (model.hasOptionalEventOccured());
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
//	        	reactionCounter++;
//	        	reactionCounterArray[reaction]++;
	        	for (TrajectoryRecorder<T> handler : trajectoryRecorders)
	        		handler.handleReactionEvent(reaction, t, x);
	        	model.checkForOptionalEvent(t, x);
	        }
		}
//		System.out.println("Total of " + evaluationCounter + " evaluations and " + reactionCounter + " reactions performed");
//		Utilities.printArray("Total reaction counts", reactionCounterArray);
//		for (int r=0; r < reactionCounterArray.length; r++)
//			reactionCounterArray[r] /= reactionCounter;
//		Utilities.printArray("Relative reaction counts", reactionCounterArray);
		integrator.clearEventHandlers();
		integrator.clearStepHandlers();
		for (int i=0; i < x1.length; i++)
			x1[i] = x[i];
    	for (TrajectoryRecorder<T> handler : trajectoryRecorders)
    		handler.setFinalState(t1, x1);
		return t;
	}

	@Override
	public void addTrajectoryRecorder(ContinuousTrajectoryRecorder<T> tr) {
		trajectoryRecorders.add(tr);
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