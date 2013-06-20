package ch.ethz.khammash.hybridstochasticsimulation.networks;


public interface ReactionNetwork {

	int getNumberOfSpecies();

	int getNumberOfReactions();

	int getStochiometry(int species, int reaction);

	int[] getStochiometries(int reaction);

	int[][] getStochiometries();

}
