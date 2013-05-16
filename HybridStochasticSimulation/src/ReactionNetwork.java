import java.util.ArrayList;
import java.util.List;


/*
 * The rateParameters k_i are of the form so that the propensity of a reaction is
 * 	k_i for constitutive reactions
 *  k_i * x1 for unary reactions
 *  k_i * x1 * x2 for binary reactions of different species
 *  1/2 * k_i * x1 * (x1 - 1) for binary reactions of the same species
 */
public class ReactionNetwork {

	protected int numOfSpecies;
	protected int numOfReactions;
	protected int[][] productionStochiometry;
	protected int[][] consumptionStochiometry;
	protected int[][] stochiometry;
	private double[] rateParameters;

	public ReactionNetwork(int numOfSpecies, int numOfReactions) {
		this.numOfSpecies = numOfSpecies;
		this.numOfReactions = numOfReactions;
		productionStochiometry = new int[numOfReactions][numOfSpecies];
		consumptionStochiometry = new int[numOfReactions][numOfSpecies];
		stochiometry = new int[numOfReactions][numOfSpecies];
		rateParameters = new double[numOfReactions];
	}

	public int getNumberOfSpecies() {
		return numOfSpecies;
	}

	public int getNumberOfReactions() {
		return numOfReactions;
	}

	public void setStochiometry(int species, int reaction, int production, int consumption) {
		productionStochiometry[reaction][species] = production;
		consumptionStochiometry[reaction][species] = consumption;
		stochiometry[reaction][species] = production - consumption;
	}

	public void setStochiometries(int[][] productionStoch, int[][] consumptionStoch) {
		if (productionStoch.length != numOfReactions)
			throw new IndexOutOfBoundsException();
		if (productionStoch.length > 0 && productionStoch[0].length != numOfSpecies)
			throw new IndexOutOfBoundsException();
		if (consumptionStoch.length != numOfReactions)
			throw new IndexOutOfBoundsException();
		if (consumptionStoch.length > 0 && consumptionStoch[0].length != numOfSpecies)
			throw new IndexOutOfBoundsException();
		for (int r=0; r < numOfReactions; r++)
			for (int s=0; s < numOfSpecies; s++) {
				productionStochiometry[r][s] = productionStoch[r][s];
				consumptionStochiometry[r][s] = consumptionStoch[r][s];
				stochiometry[r][s] = productionStoch[r][s] - consumptionStoch[r][s];
			}
	}

	public int getProductionStochiometry(int species, int reaction) {
		return productionStochiometry[reaction][species];
	}

	public int getConsumptionStochiometry(int species, int reaction) {
		return productionStochiometry[reaction][species];
	}

	public int getStochiometry(int species, int reaction) {
		return stochiometry[reaction][species];
	}

	public int[] getProductionStochiometries(int reaction) {
		return productionStochiometry[reaction].clone();
	}

	public int[] getConsumptionStochiometries(int reaction) {
		return consumptionStochiometry[reaction].clone();
	}

	public int[] getStochiometries(int reaction) {
		return stochiometry[reaction].clone();
	}

	public int[][] getProductionStochiometries() {
		return productionStochiometry.clone();
	}

	public int[][] getConsumptionStochiometries() {
		return consumptionStochiometry.clone();
	}

	public int[][] getStochiometries() {
		return stochiometry.clone();
	}

	public void setRateParameter(int reaction, double parameter) {
		this.rateParameters[reaction] = parameter;
	}

	public void setRateParameters(double[] rateParameters) {
		if (rateParameters.length != numOfReactions)
			throw new IndexOutOfBoundsException();
		for (int r=0; r < numOfReactions; r++) {
			this.rateParameters[r] = rateParameters[r];
		}
	}

	public double getRateParameter(int reaction) {
		return rateParameters[reaction];
	}

	public double[] getRateParameters() {
		return rateParameters.clone();
	}

	public List<int[]> getChoiceIndices() {
		List<int[]> result = new ArrayList<int[]>();
		for (int r=0; r < numOfReactions; r++)
			result.add(getChoiceIndices(r));
		return result;
	}
	public int[] getChoiceIndices(int reaction) {
		int s1 = -1;
		int s2 = -1;
		for (int s=0; s < numOfSpecies; s++)
			switch (consumptionStochiometry[reaction][s]) {
			case 0:
				break;
			case 1:
				if (s1 == -1)
					s1 = s;
				else if (s2 == -1)
					s2 = s;
				else
					throw new RuntimeException(
							"Only constitutive, unary and binary reactions are allowed");
				break;
			case 2:
				if (s1 == -1) {
					s1 = s;
					s2 = s;
				} else
					throw new RuntimeException(
							"Only constitutive, unary and binary reactions are allowed");
				break;
			default:
				throw new RuntimeException(
						"Only constitutive, unary and binary reactions are allowed");
			}
		int[] result = null;
		if (s1 != -1 && s2 != -1) {
			result = new int[2];
			result[0] = s1;
			result[1] = s2;
		} else if (s1 != -1) {
			result = new int[1];
			result[0] = s1;
		} else if (s2 != -1) {
			result = new int[1];
			result[0] = s2;
		} else {
			result = new int[0];
		}
		return result;
	}

}
