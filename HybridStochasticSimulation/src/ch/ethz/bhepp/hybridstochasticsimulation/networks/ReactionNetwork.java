/**
 * @author      Benjamin Hepp <benjamin.hepp@bsse.ethz.ch>
 */
package ch.ethz.bhepp.hybridstochasticsimulation.networks;

import java.util.List;
import java.util.Set;

import ch.ethz.bhepp.hybridstochasticsimulation.graphs.ReactionNetworkGraph;

/**
 * An Interface describing a chemical reaction network.
 */
public interface ReactionNetwork {

	/**
	 * @return Nspecies, the number of chemical species in the network
	 */
	int getNumberOfSpecies();

	/**
	 * @return Nreactions, the number of reactions in the network
	 */
	int getNumberOfReactions();

	/**
	 * @param s the index of a chemical species
	 * @param r the index of a reaction
	 * @return the stoichiometry for species s and reaction r
	 */
	int getStoichiometry(int s, int r);

	/**
	 * @param r the index of a reaction
	 * @return the stoichiometry vector of reaction r
	 */
	int[] getStoichiometries(int r);

	/**
	 * @return the stoichiometry matrix as a two-dimensional array (Nspecies x Nreactions)
	 */
	int[][] getStoichiometries();

	/**
	 * @param r the index of a reaction
	 * @return the set of species modified by reaction r
	 */
	Set<Integer> getModifiedSpecies(int r);

	/**
	 * @param s the index of a chemical species
	 * @return the set of reactions modifying the species s
	 */
	Set<Integer> getModifyingReactions(int s);

	/**
	 * @return the {@link ReactionNetworkGraph} corresponding to this reaction network
	 */
	ReactionNetworkGraph getGraph();

	/**
	 * @param s the index of a chemical species
	 * @return the name of the species s
	 */
	String getSpeciesLabel(int s);

	/**
	 * @param r the index of a reaction
	 * @return the name of reaction r
	 */
	String getReactionLabel(int r);

	/**
	 * @return a list with a name for each chemical species
	 */
	List<String> getSpeciesLabels();

	/**
	 * @return a list with a name for each reaction
	 */
	List<String> getReactionLabels();

	/**
	 * @return the initial condition vector of the chemical species
	 */
	double[] getInitialConditions();

}
