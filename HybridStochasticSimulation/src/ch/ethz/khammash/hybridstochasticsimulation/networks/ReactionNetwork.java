package ch.ethz.khammash.hybridstochasticsimulation.networks;

import java.util.List;


public interface ReactionNetwork {

	int getNumberOfSpecies();

	int getNumberOfReactions();

	int getStochiometry(int species, int reaction);

	int[] getStochiometries(int reaction);

	int[][] getStochiometries();

	int getProductionStochiometry(int species, int reaction);

	int getConsumptionStochiometry(int species, int reaction);

	int[] getProductionStochiometries(int reaction);

	int[] getConsumptionStochiometries(int reaction);

	int[][] getProductionStochiometries();

	int[][] getConsumptionStochiometries();

	List<Integer> getInvolvedSpecies(int reaction);

	List<Integer> getInvolvedReactions(int species);

}
