package ch.ethz.khammash.hybridstochasticsimulation.simulators;

import java.util.ArrayList;
import java.util.Collection;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.events.EventHandler;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;


public class PDMPModelSimulator extends StochasticModelSimulator{

	protected FirstOrderIntegrator integrator;
	protected Collection<PDMPStepHandler> stepHandlers;
	private double ehMaxCheckInterval;
	private Double ehConvergence;
	private Double ehConvergenceFactor;
	private int ehMaxIterationCount;

	public PDMPModelSimulator() {
		this(null, null);
	}

	public PDMPModelSimulator(RandomDataGenerator rng) {
		this(null, rng);
	}

	public PDMPModelSimulator(FirstOrderIntegrator integrator) {
		this(integrator, null);
	}

	public PDMPModelSimulator(FirstOrderIntegrator integrator, RandomDataGenerator rng) {
		super(rng);
		ehMaxCheckInterval = Double.POSITIVE_INFINITY;
		ehConvergence = Double.valueOf(1e-10);
		ehConvergenceFactor = null;
		ehMaxIterationCount = 1000;
		if (integrator == null)
			integrator = new DormandPrince853Integrator(1.0e-8, 100.0, 1.0e-10, 1.0e-10);
		else
			this.integrator = integrator;
		stepHandlers = new ArrayList<PDMPStepHandler>();
	}

	public void setEventHandlerMaxCheckInterval(double maxCheckInterval) {
		this.ehMaxCheckInterval = maxCheckInterval;
	}

	public void setEventHandlerConvergence(double convergence) {
		this.ehConvergence = Double.valueOf(convergence);
		this.ehConvergenceFactor = null;
	}

	public void setEventHandlerConvergenceFactor(double convergenceFactor) {
		this.ehConvergence = null;
		this.ehConvergenceFactor = Double.valueOf(convergenceFactor);
	}

	public void setEventHandlerMaxIterationCount(int maxIterationCount) {
		this.ehMaxIterationCount = maxIterationCount;
	}

	public double simulate(PDMPModel model, double t0, double[] x0, double t1, double[] x1) {
    	FirstOrderDifferentialEquations ode = model.getFirstOrderDifferentialEquations();
    	ReactionNetworkModel rnm = model.getReactionNetworkModel();
    	EventHandler pdmpEventHandler = model.getPDMPEventHandler();
    	model.initialize(t0, x0);
		double[] x = new double[x0.length + 2];
		for (int i=0; i < x0.length; i++)
			x[i] = x0[i];
		double[] xDot = x.clone();
		double t = t0;
		double[] propVec = new double[rnm.getPropensityDimension()];
		for (PDMPStepHandler handler : stepHandlers) {
			handler.reset();
			handler.init(t0, x, t1);
			integrator.addStepHandler(handler);
		}
		double conv = (ehConvergence != null) ? ehConvergence : ehConvergenceFactor * (t1 - t0);
		integrator.addEventHandler(pdmpEventHandler, ehMaxCheckInterval, conv, ehMaxIterationCount);
    	for (ReactionEventHandler handler : reactionHandlers)
    		handler.setInitialState(t0, x0);
    	for (EventHandler eh : model.getOptionalEventHandlers())
    		integrator.addEventHandler(eh, ehMaxCheckInterval, conv, ehMaxIterationCount);
		while (t < t1) {
			boolean hasDeterministicPart = model.hasDeterministicPart();
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
		        x[x.length - 1] = -Math.log(rng.nextUniform(0.0,  1.0));
		        do {
		        	model.handleOptionalEvent(t, x);
		        	t = integrator.integrate(ode, t, x, t1, x);
		        } while (model.getOptionalEventFlag());
			} else {
		        rnm.computePropensities(t, x, propVec);
		        for (int i=0; i < propVec.length; i++)
		        	propSum += propVec[i];
		        propensitiesComputed = true;
		        // Find next reaction time point
		        if (propSum == 0)
		        	break;
		        // -Math.log(rng.nextUniform(0.0,  1.0))
		        double tau = rng.nextExponential(1 / propSum);
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
	        double u = rng.nextUniform(0.0, 1.0);
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
	        	for (ReactionEventHandler handler : reactionHandlers)
	        		handler.handleReactionEvent(reaction, t, x);
	        	model.manualCheckOptionalEvent(t, x);
	        }
		}
		integrator.clearEventHandlers();
		integrator.clearStepHandlers();
		for (int i=0; i < x1.length; i++)
			x1[i] = x[i];
    	for (ReactionEventHandler handler : reactionHandlers)
    		handler.setFinalState(t1, x1);
		return t;
	}

	public void addStepHandler(PDMPStepHandler handler) {
		if (stepHandlers.contains(handler) == false)
			stepHandlers.add(handler);
	}
	public void removeStepHandler(PDMPStepHandler handler) {
		if (stepHandlers.contains(handler))
			stepHandlers.remove(handler);
	}
	public Collection<PDMPStepHandler> getStepHandlers() {
		return stepHandlers;
	}
	public void clearStepHandlers() {
		stepHandlers.clear();
	}

}
