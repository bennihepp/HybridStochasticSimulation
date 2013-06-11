package ch.ethz.khammash.hybridstochasticsimulation.networks;

import java.util.List;

public interface UnaryBinaryReactionNetwork extends ReactionNetwork {

	public int getProductionStochiometry(int species, int reaction);

	public int getConsumptionStochiometry(int species, int reaction);

	public int[] getProductionStochiometries(int reaction);

	public int[] getConsumptionStochiometries(int reaction);

	public int[][] getProductionStochiometries();

	public int[][] getConsumptionStochiometries();

	public double getRateParameter(int reaction);

	public double[] getRateParameters();

	public int[] getChoiceIndices(int reaction);

	public List<int[]> getChoiceIndicesList();

}
