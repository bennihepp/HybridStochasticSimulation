package ch.ethz.khammash.hybridstochasticsimulation.models;

import static com.google.common.base.Preconditions.checkElementIndex;

public class HybridReactionNetwork extends ReactionNetwork {

	public enum SpeciesType {
		DISCRETE, CONTINUOUS,
	}

	public enum ReactionType {
		NONE, STOCHASTIC, DETERMINISTIC, HYBRID,
	}

	private SpeciesType[] speciesTypes;

	public HybridReactionNetwork(int numOfSpecies, int numOfReactions, int[] continuousSpecies) {
		super(numOfSpecies, numOfReactions);
		this.setDiscreteSpecies(continuousSpecies);
	}

	public HybridReactionNetwork(ReactionNetwork net, int[] continuousSpecies) {
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
		boolean stochasticReaction = true;
		boolean deterministicReaction = true;
		for (int s=0; s < getNumberOfSpecies(); s++) {
			SpeciesType sType = getSpeciesType(s);
			int stochiometry = getStochiometry(s, r);
			if (stochiometry != 0) {
				if (sType == SpeciesType.DISCRETE)
					deterministicReaction = false;
				else if (sType == SpeciesType.CONTINUOUS)
					stochasticReaction = false;
			}
		}
		ReactionType rType = ReactionType.NONE;
		if (stochasticReaction && deterministicReaction)
			rType = ReactionType.NONE;
		else if (!stochasticReaction && !deterministicReaction)
			rType = ReactionType.HYBRID;
		else if (stochasticReaction && !deterministicReaction)
			rType = ReactionType.STOCHASTIC;
		else if (!stochasticReaction && deterministicReaction)
			rType = ReactionType.DETERMINISTIC;
		return rType;
	}

}
