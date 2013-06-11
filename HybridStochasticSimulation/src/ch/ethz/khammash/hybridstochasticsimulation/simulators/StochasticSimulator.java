package ch.ethz.khammash.hybridstochasticsimulation.simulators;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.ArrayList;
import java.util.Collection;

import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;


public class StochasticSimulator {

	protected RandomDataGenerator rdg;
	protected Collection<ReactionEventHandler> reactionHandlers;

	public StochasticSimulator() {
		this(null);
	}

	public StochasticSimulator(RandomDataGenerator rdg) {
		if (rdg == null)
			rdg = new RandomDataGenerator();
		this.rdg = rdg;
		reactionHandlers = new ArrayList<ReactionEventHandler>();
	}

	public double simulate(StochasticReactionNetworkModel model, double t0, double[] x0, double t1, double[] x1) {
		checkArgument(x0.length == x1.length, "Expected x0.length == x1.length");
		checkArgument(x0.length == model.getNumberOfSpecies(), "Expected x0.length == model.getNumberOfSpecies()");
		double[] x = new double[x0.length];
		for (int i=0; i < x0.length; i++)
			x[i] = x0[i];
		double t = t0;
		double[] propVec = new double[model.getNumberOfReactions()];
    	for (ReactionEventHandler handler : reactionHandlers)
    		handler.setInitialState(t0, x0);
//    	long reactionCounter = 0;
//    	double[] reactionCounterArray = new double[model.getPropensityDimension()];
		while (true) {
	        model.computePropensities(t, x, propVec);
	        double propSum = 0.0;
	        for (int i=0; i < propVec.length; i++)
	        	propSum += propVec[i];
	        // Find next reaction time point
	        if (propSum == 0)
	        	break;
	        double tau = rdg.nextExponential(1 / propSum);
	        t = t + tau;
	        if (t >= t1)
	        	break;
	        // Determine which reaction fired and update state
	        double u = rdg.nextUniform(0.0, 1.0);
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
	        if (reaction >= 0) {
//	        	reactionCounter++;
//	        	reactionCounterArray[reaction]++;
	        	for (ReactionEventHandler handler : reactionHandlers)
	        		handler.handleReactionEvent(reaction, t, x);
	        }
		}
//		System.out.println("Total of " + reactionCounter + " reactions performed");
//		Utilities.printArray("Total reaction counts", reactionCounterArray);
//		for (int r=0; r < reactionCounterArray.length; r++)
//			reactionCounterArray[r] /= reactionCounter;
//		Utilities.printArray("Relative reaction counts", reactionCounterArray);
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
