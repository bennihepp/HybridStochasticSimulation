/**
 * @author      Benjamin Hepp <benjamin.hepp@bsse.ethz.ch>
 */
package ch.ethz.bhepp.hybridstochasticsimulation.networks;

import java.util.List;
import java.util.Set;


/**
 * An Interface describing a chemical reaction network with mass action kinetics.
 * 
 * The interface provides methods to access the production and consumption stochiometries
 * (stoichiometry = production stoichiometry - consumption stoichiometry) and to access
 * rate parameters and details of reactants and products.
 */
public interface MassActionReactionNetwork extends ReactionNetwork {

	/**
	 * @param s the index of a chemical species
	 * @param r the index of a reaction
	 * @return the production stoichiometry for species s and reaction r
	 */
	int getProductionStoichiometry(int s, int r);

	/**
	 * @param s the index of a chemical species
	 * @param r the index of a reaction
	 * @return the consumption stoichiometry for species s and reaction r
	 */
	int getConsumptionStoichiometry(int s, int r);

	/**
	 * @param r the index of a reaction
	 * @return the production stoichiometry vector of reaction r
	 */
	int[] getProductionStoichiometries(int r);

	/**
	 * @param r the index of a reaction
	 * @return the consumption stoichiometry vector of reaction r
	 */
	int[] getConsumptionStoichiometries(int r);

	/**
	 * @return the production stoichiometry matrix as a two-dimensional array (Nspecies x Nreactions)
	 */
	int[][] getProductionStoichiometries();

	/**
	 * @return the consumption stoichiometry matrix as a two-dimensional array (Nspecies x Nreactions)
	 */
	int[][] getConsumptionStoichiometries();

	/**
	 * @param r the index of a reaction
	 * @return the set of species that are involved in reaction r, i.e. production stoichiometry != 0 or consumption stoichiometry != 0
	 */
	Set<Integer> getInvolvedSpecies(int r);

	/**
	 * @param s the index of a chemical species
	 * @return the set of reactions involving species s, i.e. production stoichiometry != 0 or consumption stoichiometry != 0
	 */
	Set<Integer> getInvolvingReactions(int s);

	/**
	 * @param r the index of a reaction
	 * @return the rate parameter of reaction r
	 */
	double getRateParameter(int r);

	/**
	 * @return the vector of rate parameters
	 */
	double[] getRateParameters();

	/**
	 * For a chemical reaction 2A + B -> ... getReactantIndices will yield [0, 0, 1] (the order is not guaranteed) if
	 * species A has index 0 and species B has index 1.
	 * 
	 * @param r the index of a reaction
	 * @return the indices of reactants in reaction r
	 */
	int[] getReactantIndices(int r);

	/**
	 * For a chemical reaction 2A + B -> ... getReactantIndices will yield [0, 0, 1] (the order is not guaranteed) if
	 * species A has index 0 and species B has index 1.
	 * 
	 * @param r the index of a reaction
	 * @return the indices of reactants in reaction r
	 */
	List<Integer> getReactantIndicesAsList(int r);

}
