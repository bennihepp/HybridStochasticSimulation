package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.io.Serializable;
import java.util.Arrays;

import ch.ethz.khammash.hybridstochasticsimulation.networks.ReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;


public class UnaryBinaryStochasticModel implements StochasticReactionNetworkModel, Serializable {

	private static final long serialVersionUID = 1L;

	final private int numberOfSpecies;
	private int[] reactionChoiceIndices1;
	private int[] reactionChoiceIndices2;
	private double[] modelRateParameters;
	private int[][] reactionStochiometries;
	private UnaryBinaryReactionNetwork network;

    public UnaryBinaryStochasticModel(UnaryBinaryReactionNetwork network) {
    	this.network = network;
    	numberOfSpecies = network.getNumberOfSpecies();
    	init(network);
    }

    final private void init(UnaryBinaryReactionNetwork net) {
    	reactionChoiceIndices1 = new int[net.getNumberOfReactions()];
    	reactionChoiceIndices2 = new int[net.getNumberOfReactions()];
    	for (int r=0; r < net.getNumberOfReactions(); r++) {
    		reactionChoiceIndices1[r] = -1;
    		reactionChoiceIndices2[r] = -1;
    	}
		modelRateParameters = net.getRateParameters();
    	reactionStochiometries = net.getStoichiometries();

    	for (int r=0; r < net.getNumberOfReactions(); r++) {
    		int[] choiceIndices = net.getChoiceIndices(r);
    		switch (choiceIndices.length) {
    		case 0:
    			break;
    		case 1:
    			reactionChoiceIndices1[r] = choiceIndices[0];
    			break;
    		case 2:
    			reactionChoiceIndices1[r] = choiceIndices[0];
    			reactionChoiceIndices2[r] = choiceIndices[1];
    			if (choiceIndices[0] == choiceIndices[1])
    				modelRateParameters[r] *= 1 / 2.0;
    			break;
			default:
				throw new RuntimeException("Only constitutive, unary and binary reactions are allowed");
    		}
    	}
    }

    @Override
    public int getNumberOfSpecies() {
    	return numberOfSpecies;
    }

	@Override
	public int getNumberOfReactions() {
		return modelRateParameters.length;
	}

	@Override
	final public double computePropensity(int reaction, double t, double[] x) {
		double p = modelRateParameters[reaction];
		int choiceIndex1 = reactionChoiceIndices1[reaction];
		int choiceIndex2 = reactionChoiceIndices2[reaction];
		if (choiceIndex1 != -1) {
			p *= x[choiceIndex1];
			if (choiceIndex2 != -1) {
				double q = x[choiceIndex2];
				if (choiceIndex1 == choiceIndex2)
					q -= 1;
				p *= q;
			}
		}
		return p;
	}

	@Override
	public void computePropensities(double t, double[] x, double[] propensities) {
		// We don't check the length of x and propensities for performance reasons
		Arrays.fill(propensities, 0, getNumberOfReactions(), 0.0);
    	for (int r=0; r < getNumberOfReactions(); r++)
    		propensities[r] = computePropensity(r, t, x);
	}

	@Override
	public double computePropensitiesAndSum(double t, double[] x, double[] propensities) {
		double propSum = 0.0;
		// We don't check the length of x and propensities for performance reasons
		Arrays.fill(propensities, 0, getNumberOfReactions(), 0.0);
    	for (int r=0; r < getNumberOfReactions(); r++) {
    		double propensity = computePropensity(r, t, x); 
    		propensities[r] = propensity;
    		propSum += propensity;
    	}
    	return propSum;
	}

	@Override
	public void changeState(int reaction, double t, double[] x) {
		// We don't check the length of x and propensities for performance reasons
		int[] stochiometry = reactionStochiometries[reaction];
    	for (int s=0; s < stochiometry.length; s++) {
    		x[s] += stochiometry[s];
    	}
	}

	@Override
	public boolean isTimeIndependent() {
		return true;
	}

	@Override
	public ReactionNetwork getNetwork() {
		return network;
	}

}
