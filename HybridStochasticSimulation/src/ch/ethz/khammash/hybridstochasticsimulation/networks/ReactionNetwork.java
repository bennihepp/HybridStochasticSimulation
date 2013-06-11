package ch.ethz.khammash.hybridstochasticsimulation.networks;


public interface ReactionNetwork {

	public int getNumberOfSpecies();

	public int getNumberOfReactions();

	public int getStochiometry(int species, int reaction);

	public int[] getStochiometries(int reaction);

	public int[][] getStochiometries();

}
