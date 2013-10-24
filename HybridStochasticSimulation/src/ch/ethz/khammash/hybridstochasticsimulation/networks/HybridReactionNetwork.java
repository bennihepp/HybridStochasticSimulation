package ch.ethz.khammash.hybridstochasticsimulation.networks;

import static com.google.common.base.Preconditions.checkElementIndex;

import java.io.Serializable;
import java.util.Arrays;

public class HybridReactionNetwork extends DefaultUnaryBinaryReactionNetwork implements Serializable {

	private static final long serialVersionUID = 1L;


	public enum ReactionType {
		STOCHASTIC, DETERMINISTIC,
	}

	private ReactionType[] reactionTypes;

	public HybridReactionNetwork(int numOfSpecies, int numOfReactions) {
		this(numOfSpecies, numOfReactions, new int[0]);
	}

	public HybridReactionNetwork(int numOfSpecies, int numOfReactions, int[] deterministicReactions) {
		super(numOfSpecies, numOfReactions);
		initReactionTypes();
		this.setDeterministicReactions(deterministicReactions);
	}

	public HybridReactionNetwork(UnaryBinaryReactionNetwork net) {
		this(net, new int[0]);
	}

	public HybridReactionNetwork(UnaryBinaryReactionNetwork net, int[] deterministicReactions) {
		super(net.getNumberOfSpecies(), net.getNumberOfReactions());
		setStoichiometries(net.getProductionStoichiometries(), net.getConsumptionStoichiometries());
		setRateParameters(net.getRateParameters());
		initReactionTypes();
		this.setDeterministicReactions(deterministicReactions);
	}

	private final void initReactionTypes() {
		reactionTypes = new ReactionType[getNumberOfReactions()];
		Arrays.fill(reactionTypes, ReactionType.STOCHASTIC);
	}

	public void setDeterministicReactions(int[] deterministicReactions) {
		for (int r : deterministicReactions) {
			checkElementIndex(r, getNumberOfReactions(), "Expected 0 <= deterministicReactions[i] < getNumberOfReactions()");
			reactionTypes[r] = ReactionType.DETERMINISTIC;
		}
	}

	public void setDeterministicReaction(int reaction, int deterministicReactions) {
		checkElementIndex(reaction, getNumberOfReactions(), "Expected 0 <= reaction < getNumberOfReactions()");
		reactionTypes[reaction] = ReactionType.DETERMINISTIC;
	}


	public ReactionType getReactionType(int r) {
		checkElementIndex(r, getNumberOfReactions(), "Expected 0 <= s < getNumberOfReactions()");
		return reactionTypes[r];
	}

}
