package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.Arrays;

import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;


public class UnaryBinaryStochasticModel implements StochasticReactionNetworkModel {

	final int numberOfSpecies;
	protected int[] reactionChoiceIndex1;
	protected int[] reactionChoiceIndex2;
	protected double[] rateParameters;
	protected int[][] reactionStochiometries;

    public UnaryBinaryStochasticModel(UnaryBinaryReactionNetwork net) {
    	numberOfSpecies = net.getNumberOfSpecies();
    	init(net);
    }

    final private void init(UnaryBinaryReactionNetwork net) {
    	reactionChoiceIndex1 = new int[net.getNumberOfReactions()];
    	reactionChoiceIndex2 = new int[net.getNumberOfReactions()];
    	for (int r=0; r < net.getNumberOfReactions(); r++) {
    		reactionChoiceIndex1[r] = -1;
    		reactionChoiceIndex2[r] = -1;
    	}
    	rateParameters = net.getRateParameters();
    	reactionStochiometries = net.getStochiometries();

    	for (int r=0; r < net.getNumberOfReactions(); r++) {
    		int[] choiceIndices = net.getChoiceIndices(r);
    		switch (choiceIndices.length) {
    		case 0:
    			break;
    		case 1:
    			reactionChoiceIndex1[r] = choiceIndices[0];
    			break;
    		case 2:
    			reactionChoiceIndex1[r] = choiceIndices[0];
    			reactionChoiceIndex2[r] = choiceIndices[1];
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
		return rateParameters.length;
	}

	@Override
	final public double computePropensity(int reaction, double t, double[] x) {
		// We don't check the length of x and propensities for performance reasons
		double p = rateParameters[reaction];
		int choiceIndex1 = reactionChoiceIndex1[reaction];
		int choiceIndex2 = reactionChoiceIndex2[reaction];
		if (choiceIndex2 != -1) {
			if (choiceIndex1 == choiceIndex2)
				// Binary reaction of the same species
				p *= (1/2.0) * x[choiceIndex1] * (x[choiceIndex1] - 1);
			else
				// Binary reaction of two different species
				p *= x[choiceIndex1] * x[choiceIndex2];
		} else if (choiceIndex1 != -1) {
			// Unary reaction
			p *= x[choiceIndex1];
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
	public void updateState(int reaction, double t, double[] x) {
		// We don't check the length of x and propensities for performance reasons
		int[] stochiometry = reactionStochiometries[reaction];
    	for (int i=0; i < stochiometry.length; i++) {
    		x[i] += stochiometry[i];
    	}
	}

	@Override
	public boolean isTimeIndependent() {
		return true;
	}

}
