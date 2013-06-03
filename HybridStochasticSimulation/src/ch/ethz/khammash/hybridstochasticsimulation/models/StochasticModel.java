package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.List;

import ch.ethz.khammash.hybridstochasticsimulation.networks.ReactionNetwork;


// TODO: Group computation of choices for each reaction

public class StochasticModel implements ReactionNetworkModel {

	final int numberOfStates;
	protected int[] reactionChoiceIndex1;
	protected int[] reactionChoiceIndex2;
	protected double[] rateParameters;
	protected int[][] reactionStochiometries;

    public StochasticModel(ReactionNetwork net) {
    	numberOfStates = net.getNumberOfSpecies();
    	init(net);
    }

    protected void init(ReactionNetwork net) {
    	List<int[]> choiceIndicesList = net.getChoiceIndices();

    	reactionChoiceIndex1 = new int[net.getNumberOfReactions()];
    	reactionChoiceIndex2 = new int[net.getNumberOfReactions()];
    	for (int r=0; r < net.getNumberOfReactions(); r++) {
    		reactionChoiceIndex1[r] = -1;
    		reactionChoiceIndex2[r] = -1;
    	}
    	rateParameters = net.getRateParameters();
    	reactionStochiometries = net.getStochiometries();

    	for (int r=0; r < net.getNumberOfReactions(); r++) {
    		int[] choiceIndices = choiceIndicesList.get(r);
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
    public int getNumberOfStates() {
    	return numberOfStates;
    }

	@Override
	public int getPropensityDimension() {
		return rateParameters.length;
	}

	@Override
	public void computePropensities(double t, double[] x, double[] propensities) {
		// We don't check the length of x and propensities for performance reasons
    	for (int i=0; i < rateParameters.length; i++) {
    		double p = rateParameters[i];
    		int choiceIndex1 = reactionChoiceIndex1[i];
    		int choiceIndex2 = reactionChoiceIndex2[i];
    		if (choiceIndex2 != -1) {
    			if (choiceIndex1 == choiceIndex2)
    				// binary reaction of the same species
    				p *= (1/2.0) * x[choiceIndex1] * (x[choiceIndex1] - 1);
    			else
    				// binary reaction of two different species
    				p *= x[choiceIndex1] * x[choiceIndex2];
    		} else if (choiceIndex1 != -1) {
    			// Unary reaction
    			p *= x[choiceIndex1];
    		}
    		propensities[i] = p;
    	}
	}

	@Override
	public void updateState(int reaction, double t, double[] x) {
		// We don't check the length of x and propensities for performance reasons
		int[] stochiometry = reactionStochiometries[reaction];
    	for (int i=0; i < stochiometry.length; i++) {
    		x[i] += stochiometry[i];
    	}
	}

}
