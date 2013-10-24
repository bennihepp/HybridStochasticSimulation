package ch.ethz.khammash.hybridstochasticsimulation.networks;

import java.util.List;

public interface UnaryBinaryReactionNetwork extends ReactionNetwork {

	int getProductionStoichiometry(int species, int reaction);

	int getConsumptionStoichiometry(int species, int reaction);

	int[] getProductionStoichiometries(int reaction);

	int[] getConsumptionStoichiometries(int reaction);

	int[][] getProductionStoichiometries();

	int[][] getConsumptionStoichiometries();

	double getRateParameter(int reaction);

	double[] getRateParameters();

	int[] getChoiceIndices(int reaction);

	List<int[]> getChoiceIndicesList();

}
