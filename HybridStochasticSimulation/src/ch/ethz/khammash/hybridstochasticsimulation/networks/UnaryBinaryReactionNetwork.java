package ch.ethz.khammash.hybridstochasticsimulation.networks;

import java.util.List;

public interface UnaryBinaryReactionNetwork extends ReactionNetwork {

	int getProductionStochiometry(int species, int reaction);

	int getConsumptionStochiometry(int species, int reaction);

	int[] getProductionStochiometries(int reaction);

	int[] getConsumptionStochiometries(int reaction);

	int[][] getProductionStochiometries();

	int[][] getConsumptionStochiometries();

	double getRateParameter(int reaction);

	double[] getRateParameters();

	int[] getChoiceIndices(int reaction);

	List<int[]> getChoiceIndicesList();

}
