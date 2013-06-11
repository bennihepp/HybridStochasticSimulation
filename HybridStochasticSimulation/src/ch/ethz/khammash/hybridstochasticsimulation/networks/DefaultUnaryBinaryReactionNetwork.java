package ch.ethz.khammash.hybridstochasticsimulation.networks;

import static com.google.common.base.Preconditions.*;

import java.util.ArrayList;
import java.util.List;


/*
 * rateParameters k_i are of the form so that the propensity of a reaction is
 * 	k_i for constitutive reactions
 *  k_i * x1 for unary reactions
 *  k_i * x1 * x2 for binary reactions of different species
 *  1/2 * k_i * x1 * (x1 - 1) for binary reactions of the same species
 */
public class DefaultUnaryBinaryReactionNetwork implements UnaryBinaryReactionNetwork {

	private int numOfSpecies;
	private int numOfReactions;
	private int[][] productionStochiometry;
	private int[][] consumptionStochiometry;
	private int[][] stochiometry;
	private double[] rateParameters;
	private List<int[]> choiceIndicesList;

	public DefaultUnaryBinaryReactionNetwork(int numOfSpecies, int numOfReactions) {
		checkArgument(numOfSpecies > 0, "Expected numOfSpecies > 0");
		checkArgument(numOfReactions > 0, "Expected numOfReactions > 0");
		this.numOfSpecies = numOfSpecies;
		this.numOfReactions = numOfReactions;
		productionStochiometry = new int[numOfReactions][numOfSpecies];
		consumptionStochiometry = new int[numOfReactions][numOfSpecies];
		stochiometry = new int[numOfReactions][numOfSpecies];
		rateParameters = new double[numOfReactions];
	}

	@Override
	final public int getNumberOfSpecies() {
		return numOfSpecies;
	}

	@Override
	final public int getNumberOfReactions() {
		return numOfReactions;
	}

	final protected void invalidate() {
		choiceIndicesList = null;
	}

	public void setStochiometry(int species, int reaction, int production, int consumption) {
		checkElementIndex(species, getNumberOfSpecies(), "Expected 0<=species<getNumberOfSpecies()");
		checkElementIndex(reaction, getNumberOfReactions(), "Expected 0<=species<getNumberOfSpecies()");
		checkArgument(production >= 0, "Expected production >= 0");
		checkArgument(consumption >= 0, "Expected consumption >= 0");
		productionStochiometry[reaction][species] = production;
		consumptionStochiometry[reaction][species] = consumption;
		stochiometry[reaction][species] = production - consumption;
		invalidate();
	}

	public void setStochiometries(int[][] productionStoch, int[][] consumptionStoch) {
		checkArgument(productionStoch.length == getNumberOfReactions());
		checkArgument(consumptionStoch.length == getNumberOfReactions());
		checkArgument(productionStoch[0].length == getNumberOfSpecies());
		checkArgument(consumptionStoch[0].length == getNumberOfSpecies());
		for (int r=0; r < getNumberOfReactions(); r++) {
			System.arraycopy(productionStoch[r], 0, productionStochiometry[r], 0, numOfSpecies);
			System.arraycopy(consumptionStoch[r], 0, consumptionStochiometry[r], 0, numOfSpecies);
			for (int s=0; s < getNumberOfSpecies(); s++)
				stochiometry[r][s] = productionStoch[r][s] - consumptionStoch[r][s];
		}
		invalidate();
	}

	@Override
	final public int getProductionStochiometry(int species, int reaction) {
		checkElementIndex(species, getNumberOfSpecies());
		checkElementIndex(reaction, getNumberOfReactions());
		return productionStochiometry[reaction][species];
	}

	@Override
	final public int getConsumptionStochiometry(int species, int reaction) {
		checkElementIndex(species, getNumberOfSpecies());
		checkElementIndex(reaction, getNumberOfReactions());
		return consumptionStochiometry[reaction][species];
	}

	@Override
	final public int getStochiometry(int species, int reaction) {
		checkElementIndex(species, getNumberOfSpecies());
		checkElementIndex(reaction, getNumberOfReactions());
		return stochiometry[reaction][species];
	}

	@Override
	public int[] getProductionStochiometries(int reaction) {
		checkElementIndex(reaction, getNumberOfReactions());
		return productionStochiometry[reaction].clone();
	}

	@Override
	public int[] getConsumptionStochiometries(int reaction) {
		checkElementIndex(reaction, getNumberOfReactions());
		return consumptionStochiometry[reaction].clone();
	}

	@Override
	public int[] getStochiometries(int reaction) {
		checkElementIndex(reaction, getNumberOfReactions());
		return stochiometry[reaction].clone();
	}

	@Override
	public int[][] getProductionStochiometries() {
		return productionStochiometry.clone();
	}

	@Override
	public int[][] getConsumptionStochiometries() {
		return consumptionStochiometry.clone();
	}

	@Override
	public int[][] getStochiometries() {
		return stochiometry.clone();
	}

	public void setRateParameter(int reaction, double parameter) {
		checkElementIndex(reaction, getNumberOfReactions());
		this.rateParameters[reaction] = parameter;
	}

	public void setRateParameters(double[] rateParameters) {
		System.arraycopy(rateParameters, 0, this.rateParameters, 0, getNumberOfReactions());
		for (int r=0; r < numOfReactions; r++)
			this.rateParameters[r] = rateParameters[r];
	}

	@Override
	final public double getRateParameter(int reaction) {
		checkElementIndex(reaction, getNumberOfReactions());
		return rateParameters[reaction];
	}

	@Override
	public double[] getRateParameters() {
		return rateParameters.clone();
	}

	@Override
	public int[] getChoiceIndices(int reaction) {
		checkElementIndex(reaction, getNumberOfReactions());
		if (choiceIndicesList == null)
			computeChoiceIndices();
		return choiceIndicesList.get(reaction);
	}

	@Override
	public List<int[]> getChoiceIndicesList() {
		if (choiceIndicesList == null) {
			computeChoiceIndices();
		}
		return choiceIndicesList;
	}

	private void computeChoiceIndices() {
		choiceIndicesList = new ArrayList<int[]>();
		for (int r=0; r < getNumberOfReactions(); r++) {
			int s1 = -1;
			int s2 = -1;
			for (int s=0; s < getNumberOfSpecies(); s++)
				switch (consumptionStochiometry[r][s]) {
				case 0:
					break;
				case 1:
					if (s1 == -1)
						s1 = s;
					else if (s2 == -1)
						s2 = s;
					else
						throw new RuntimeException("Only constitutive, unary and binary reactions are allowed");
					break;
				case 2:
					if (s1 == -1) {
						s1 = s;
						s2 = s;
					} else
						throw new RuntimeException("Only constitutive, unary and binary reactions are allowed");
					break;
				default:
					throw new RuntimeException("Only constitutive, unary and binary reactions are allowed");
				}
			int[] choiceIndices = null;
			if (s1 != -1 && s2 != -1) {
				choiceIndices = new int[2];
				choiceIndices[0] = s1;
				choiceIndices[1] = s2;
			} else if (s1 != -1) {
				choiceIndices = new int[1];
				choiceIndices[0] = s1;
			} else if (s2 != -1) {
				choiceIndices = new int[1];
				choiceIndices[0] = s2;
			} else {
				choiceIndices = new int[0];
			}
			choiceIndicesList.add(choiceIndices);
		}
	}

}
