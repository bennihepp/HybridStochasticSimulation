package ch.ethz.khammash.hybridstochasticsimulation.networks;

import static com.google.common.base.Preconditions.checkElementIndex;

public class HybridReactionNetwork extends DefaultUnaryBinaryReactionNetwork {

	public enum SpeciesType {
		DISCRETE, CONTINUOUS,
	}

	public enum ReactionType {
		NONE, STOCHASTIC, DETERMINISTIC,
	}

	private SpeciesType[] speciesTypes;

	public HybridReactionNetwork(int numOfSpecies, int numOfReactions, int[] continuousSpecies) {
		super(numOfSpecies, numOfReactions);
		this.setDiscreteSpecies(continuousSpecies);
	}

	public HybridReactionNetwork(DefaultUnaryBinaryReactionNetwork net, int[] continuousSpecies) {
		super(net.getNumberOfSpecies(), net.getNumberOfReactions());
		this.setDiscreteSpecies(continuousSpecies);
		setStochiometries(net.getProductionStochiometries(), net.getConsumptionStochiometries());
		setRateParameters(net.getRateParameters());
	}

	private void setDiscreteSpecies(int[] continuousSpecies) {
		speciesTypes = new SpeciesType[getNumberOfSpecies()];
		for (int s=0; s < getNumberOfSpecies(); s++) {
			speciesTypes[s] = SpeciesType.DISCRETE;
		}
		for (int i=0; i < continuousSpecies.length; i++) {
			int s = continuousSpecies[i];
			checkElementIndex(s, getNumberOfSpecies(), "Expected 0 <= discreteSpecies[i] < getNumberOfSpecies()");
			speciesTypes[s] = SpeciesType.CONTINUOUS;
		}
	}

	public SpeciesType getSpeciesType(int s) {
		checkElementIndex(s, getNumberOfSpecies(), "Expected 0 <= s < getNumberOfSpecies()");
		return speciesTypes[s];
	}

	public ReactionType getReactionType(int r) {
		boolean anyStochasticReaction = false;
		boolean anyDeterministicReaction = false;
		for (int s=0; s < getNumberOfSpecies(); s++) {
			SpeciesType sType = getSpeciesType(s);
			int stochiometry = getStochiometry(s, r);
			if (stochiometry != 0) {
				if (sType == SpeciesType.DISCRETE)
					anyStochasticReaction = true;
				else if (sType == SpeciesType.CONTINUOUS)
					anyDeterministicReaction = true;
			}
		}
		if (anyStochasticReaction)
			return ReactionType.STOCHASTIC;
		else if (anyDeterministicReaction)
			return ReactionType.DETERMINISTIC;
		else
			return ReactionType.NONE;
	}

}
