package ch.ethz.khammash.hybridstochasticsimulation.simulators;

import static com.google.common.base.Preconditions.*;

import java.util.ArrayList;
import java.util.Collection;

import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticModel;


public class StochasticModelSimulator {

	protected RandomDataGenerator rng;
	protected Collection<ReactionEventHandler> reactionHandlers;

	public StochasticModelSimulator() {
		this(null);
	}

	public StochasticModelSimulator(RandomDataGenerator rng) {
		if (rng == null)
			rng = new RandomDataGenerator();
		this.rng = rng;
		reactionHandlers = new ArrayList<ReactionEventHandler>();
	}

	public double simulate(StochasticModel model, double t0, double[] x0, double t1, double[] x1) {
		checkArgument(x0.length == x1.length, "Expected x0.length == x1.length");
		checkArgument(x0.length == model.getNumberOfStates(), "Expected x0.length == model.getNumberOfSpecies()");
		double[] x = new double[x0.length];
		for (int i=0; i < x0.length; i++)
			x[i] = x0[i];
		double t = t0;
		double[] propVec = new double[model.getPropensityDimension()];
    	for (ReactionEventHandler handler : reactionHandlers)
    		handler.setInitialState(t0, x0);
		while (true) {
	        model.computePropensities(t, x, propVec);
	        double propSum = 0.0;
	        for (int i=0; i < propVec.length; i++)
	        	propSum += propVec[i];
	        // Find next reaction time point
	        if (propSum == 0)
	        	break;
	        // -Math.log(rng.nextUniform(0.0,  1.0))
	        double tau = rng.nextExponential(1 / propSum);
	        t = t + tau;
	        if (t >= t1)
	        	break;
	        // Determine which reaction fired and update state
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
	        	for (ReactionEventHandler handler : reactionHandlers)
	        		handler.handleReactionEvent(reaction, t, x);
		}
		for (int i=0; i < x1.length; i++)
			x1[i] = x[i];
    	for (ReactionEventHandler handler : reactionHandlers)
    		handler.setFinalState(t1, x1);
		return t;
	}

	public void addReactionHandler(ReactionEventHandler handler) {
		if (reactionHandlers.contains(handler) == false)
			reactionHandlers.add(handler);
	}
	public void removeReactionHandler(ReactionEventHandler handler) {
		if (reactionHandlers.contains(handler))
			reactionHandlers.remove(handler);
	}
	public Collection<ReactionEventHandler> getReactionHandlers() {
		return reactionHandlers;
	}
	public void clearReactionHandlers() {
		reactionHandlers.clear();
	}

}
