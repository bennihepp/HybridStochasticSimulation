package ch.ethz.khammash.hybridstochasticsimulation.networks;

import java.util.List;

import ch.ethz.khammash.hybridstochasticsimulation.graphs.ReactionNetworkGraph;


public interface ReactionNetwork {

	int getNumberOfSpecies();

	int getNumberOfReactions();

	int getStoichiometry(int species, int reaction);

	int[] getStoichiometries(int reaction);

	int[][] getStoichiometries();

	int getProductionStoichiometry(int species, int reaction);

	int getConsumptionStoichiometry(int species, int reaction);

	int[] getProductionStoichiometries(int reaction);

	int[] getConsumptionStoichiometries(int reaction);

	int[][] getProductionStoichiometries();

	int[][] getConsumptionStoichiometries();

	List<Integer> getInvolvedSpecies(int reaction);

	List<Integer> getInvolvedReactions(int species);

	ReactionNetworkGraph getGraph();

	String getSpeciesLabel(int species);

	String getReactionLabel(int reaction);

	List<String> getSpeciesLabels();

	List<String> getReactionLabels();

	double[] getInitialConditions();

}
