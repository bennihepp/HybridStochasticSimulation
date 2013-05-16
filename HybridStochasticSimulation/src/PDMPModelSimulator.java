import java.util.ArrayList;
import java.util.Collection;

import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;
import org.apache.commons.math3.random.RandomDataGenerator;


public class PDMPModelSimulator extends StochasticModelSimulator{

	protected FirstOrderIntegrator integrator;
	protected Collection<PDMPStepHandler> stepHandlers;

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
		if (integrator == null)
			integrator = new DormandPrince853Integrator(1.0e-8, 100.0, 1.0e-10, 1.0e-10);
		if (rng == null)
			rng = new RandomDataGenerator();
		this.integrator = integrator;
		stepHandlers = new ArrayList<PDMPStepHandler>();
	}

	public double simulate(PDMPModel model, double t0, double[] x0, double t1, double[] x1) {
		double[] x = new double[x0.length + 2];
		for (int i=0; i < x0.length; i++)
			x[i] = x0[i];
		double t = t0;
		double[] propVec = new double[model.getPropensityDimension()];
		for (PDMPStepHandler handler : stepHandlers) {
			handler.reset();
			integrator.addStepHandler(handler);
		}
		integrator.addEventHandler(model, Double.POSITIVE_INFINITY, 1e-8*(t1-t0), 1000);
    	for (ReactionHandler handler : reactionHandlers)
    		handler.setInitialState(t0, x0);
		while (t < t1) {
			// Evolve ODE until next stochastic reaction fires
	        x[x.length - 1] = -Math.log(rng.nextUniform(0.0,  1.0));
	        x[x.length - 2] = 0.0;
	        t = integrator.integrate(model, t, x, t1, x);
	        // Stop if we reached the end-timepoint
	        if (t >= t1)
	        	break;

	        // Determine which reaction fired and update state
	        model.computePropensities(t, x, propVec);
	        double propSum = 0.0;
	        for (int i=0; i < propVec.length; i++)
	        	propSum += propVec[i];
	        double u = rng.nextUniform(0.0, 1.0);
	        double w = 0.0;
	        int reaction = -1;
	        for (int l=0; l < propVec.length; l++) {
	        	w = w + propVec[l] / propSum;
	        	if (u < w) {
	        		reaction = l;
	        		model.updateState(reaction, t, x);
	        		break;
	        	}
	        }
	        if (reaction >= 0)
	        	for (ReactionHandler handler : reactionHandlers)
	        		handler.handleReaction(reaction, t, x);
		}
		integrator.clearEventHandlers();
		integrator.clearStepHandlers();
		for (int i=0; i < x1.length; i++)
			x1[i] = x[i];
    	for (ReactionHandler handler : reactionHandlers)
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
