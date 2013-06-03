package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.LinkedList;
import java.util.List;

import ch.ethz.khammash.hybridstochasticsimulation.networks.HybridReactionNetwork;

// TODO: Group computation of choices for each reaction

public class HybridReactionNetworkModel implements HybridModel {

	protected int dimension;

//	protected ArrayList<int[]> stochasticReactionChoiceIndices;
//	protected double[] stochasticRateParameters;
//	protected int[][] stochasticReactionStochiometries;
//
//	protected ArrayList<int[]> deterministicReactionChoiceIndices;
//	protected double[] deterministicRateParameters;
//	protected int[][] deterministicReactionStochiometries;
//
//	protected ArrayList<int[]> hybridReactionChoiceIndices;
//	protected double[] hybridRateParameters;
//	protected int[][] hybridReactionStochiometries;

	List<int[]> choiceIndicesList;

	private int[] stochasticReactionIndices;
	private int[] deterministicReactionIndices;

	private double[] rateParameters;
	private int[] reactionChoiceIndices1;
	private int[] reactionChoiceIndices2;
	private double[][] reactionStochiometries;

	public HybridReactionNetworkModel(HybridReactionNetwork hrn) {
		choiceIndicesList = hrn.getChoiceIndices();
		dimension = hrn.getNumberOfSpecies();
		rateParameters = new double[hrn.getNumberOfReactions()];
		reactionChoiceIndices1 = new int[hrn.getNumberOfReactions()];
		reactionChoiceIndices2 = new int[hrn.getNumberOfReactions()];
		reactionStochiometries = new double[hrn.getNumberOfReactions()][hrn.getNumberOfSpecies()];

		init(hrn);
	}

//	@Override
//	public boolean hasDeterministicPart() {
//		return deterministicRateParameters.length > 0;
//	}

	@Override
	public boolean hasDeterministicPart() {
		return deterministicReactionIndices.length > 0;
	}

	private void init(HybridReactionNetwork hrn) {

		LinkedList<Integer> stochasticReactionIndicesList = new LinkedList<Integer>();
		LinkedList<Integer> deterministicReactionIndicesList = new LinkedList<Integer>();

		for (int r = 0; r < hrn.getNumberOfReactions(); r++) {
			int[] choiceIndices = choiceIndicesList.get(r);
			switch (choiceIndices.length) {
			case 0:
				reactionChoiceIndices1[r] = -1;
				reactionChoiceIndices2[r] = -1;
				break;
			case 1:
				reactionChoiceIndices1[r] = choiceIndices[0];
				reactionChoiceIndices2[r] = -1;
				break;
			case 2:
				reactionChoiceIndices1[r] = choiceIndices[0];
				reactionChoiceIndices2[r] = choiceIndices[1];
				break;
			}
			rateParameters[r] = hrn.getRateParameter(r);
			for (int s = 0; s < hrn.getNumberOfSpecies(); s++) {
				reactionStochiometries[r][s] = hrn.getStochiometry(s, r);
			}
			HybridReactionNetwork.ReactionType rt = hrn.getReactionType(r);
			switch (rt) {
			case STOCHASTIC:
				stochasticReactionIndicesList.add(Integer.valueOf(r));
				break;
			case DETERMINISTIC:
				deterministicReactionIndicesList.add(Integer.valueOf(r));
				break;
			case HYBRID:
				stochasticReactionIndicesList.add(Integer.valueOf(r));
				break;
			case NONE:
				break;
			}
		}

		stochasticReactionIndices = new int[stochasticReactionIndicesList.size()];
		deterministicReactionIndices = new int[deterministicReactionIndicesList.size()];
		for (int i=0; i < stochasticReactionIndicesList.size(); i++)
			stochasticReactionIndices[i] = stochasticReactionIndicesList.get(i).intValue();
		for (int i=0; i < deterministicReactionIndicesList.size(); i++)
			deterministicReactionIndices[i] = deterministicReactionIndicesList.get(i).intValue();

	}

//	private void init(HybridReactionNetwork hrn) {
//		int numOfStochasticReactions = 0;
//		int numOfDeterministicReactions = 0;
//		int numOfHybridReactions = 0;
//
//		for (int r = 0; r < hrn.getNumberOfReactions(); r++) {
//			HybridReactionNetwork.ReactionType rType = hrn.getReactionType(r);
//			switch (rType) {
//			case STOCHASTIC:
//				numOfStochasticReactions++;
//				break;
//			case DETERMINISTIC:
//				numOfDeterministicReactions++;
//				break;
//			case HYBRID:
//				numOfHybridReactions++;
//				break;
//			case NONE:
//				break;
//			}
//		}
//
//		stochasticReactionChoiceIndices = new ArrayList<int[]>(numOfStochasticReactions);
//		stochasticRateParameters = new double[numOfStochasticReactions];
//		stochasticReactionStochiometries = new int[numOfStochasticReactions][hrn.getNumberOfSpecies()];
//		deterministicReactionChoiceIndices = new ArrayList<int[]>(numOfDeterministicReactions);
//		deterministicRateParameters = new double[numOfDeterministicReactions];
//		deterministicReactionStochiometries = new int[numOfDeterministicReactions][hrn.getNumberOfSpecies()];
//		hybridReactionChoiceIndices = new ArrayList<int[]>(numOfHybridReactions);
//		hybridRateParameters = new double[numOfHybridReactions];
//		hybridReactionStochiometries = new int[numOfHybridReactions][hrn.getNumberOfSpecies()];
//
//		int is = 0;
//		int id = 0;
//		int ih = 0;
//
//		for (int r = 0; r < hrn.getNumberOfReactions(); r++) {
//			HybridReactionNetwork.ReactionType rType = hrn.getReactionType(r);
//			int[] choiceIndices = hrn.getChoiceIndices(r);
//			switch (rType) {
//			case STOCHASTIC:
//				stochasticReactionChoiceIndices.add(choiceIndices);
//				stochasticRateParameters[is] = hrn.getRateParameter(r);
//				stochasticReactionStochiometries[is] = hrn.getStochiometries(r);
//				is++;
//				break;
//			case DETERMINISTIC:
//				deterministicReactionChoiceIndices.add(choiceIndices);
//				deterministicRateParameters[id] = hrn.getRateParameter(r);
//				deterministicReactionStochiometries[id] = hrn.getStochiometries(r);
//				id++;
//				break;
//			case HYBRID:
//				hybridReactionChoiceIndices.add(choiceIndices);
//				hybridRateParameters[ih] = hrn.getRateParameter(r);
//				for (int s=0; s < getNumberOfStates(); s++) {
//					HybridReactionNetwork.SpeciesType sType = hrn.getSpeciesType(s);
//					switch (sType) {
//					case CONTINUOUS:
//						hybridReactionStochiometries[ih][s] = 0;
//						break;
//					case DISCRETE:
//						hybridReactionStochiometries[ih][s] = hrn.getStochiometry(s, r);
//						break;
//					}
//				}
//				ih++;
//				break;
//			case NONE:
//				break;
//			}
//		}
//	}

	@Override
	public int getDimension() {
		return dimension;
	}

	@Override
	public void computeDerivatives(double t, double[] x, double[] xDot) {
		// We don't check the length of x and propensities for performance reasons
		for (int s = 0; s < xDot.length; s++)
			xDot[s] = 0;
		for (int r : deterministicReactionIndices) {
			double v = rateParameters[r];
			int choiceIndex1 = reactionChoiceIndices1[r];
			int choiceIndex2 = reactionChoiceIndices2[r];
			if (choiceIndex1 != -1) {
				if (choiceIndex2 != -1) {
					if (choiceIndex1 == choiceIndex2)
						v *= (1 / 2.0) * x[choiceIndex1] * x[choiceIndex1];
					else
						v *= x[choiceIndex1] * x[choiceIndex2];
				} else
					v *= x[choiceIndex1];
			}
			for (int s = 0; s < reactionStochiometries[r].length; s++) {
				double stochiometry = reactionStochiometries[r][s];
				if (stochiometry != 0)
					xDot[s] += v * stochiometry;
			}
		}
	}

//	@Override
//	public void computeDerivatives(double t, double[] x, double[] xDot) {
//		// We don't check the length of x and xDot for performance reasons
//		for (int s = 0; s < xDot.length; s++)
//			xDot[s] = 0;
//		for (int i = 0; i < deterministicRateParameters.length; i++) {
//			double p = deterministicRateParameters[i];
//			int[] choiceIndices = deterministicReactionChoiceIndices.get(i);
//			if (choiceIndices.length == 1)
//				p *= x[choiceIndices[0]];
//			else if (choiceIndices.length == 2)
//				if (choiceIndices[0] == choiceIndices[1])
//					p *= (1 / 2.0) * x[choiceIndices[0]] * x[choiceIndices[1]];
//				else
//					p *= x[choiceIndices[0]] * x[choiceIndices[1]];
//			int[] stochiometry = deterministicReactionStochiometries[i];
//			for (int s=0; s < stochiometry.length; s++)
//				xDot[s] += p * stochiometry[s];
//		}
//	}

	@Override
	public int getPropensityDimension() {
		return rateParameters.length;
	}

//	@Override
//	public int getPropensityDimension() {
//		return stochasticRateParameters.length + hybridRateParameters.length;
//	}

	@Override
	public void computePropensities(double t, double[] x, double[] propensities) {
		// We don't check the length of x and propensities for performance reasons
		for (int r = 0; r < propensities.length; r++)
			propensities[r] = 0;
		for (int r : stochasticReactionIndices) {
			double p = rateParameters[r];
			int choiceIndex1 = reactionChoiceIndices1[r];
			int choiceIndex2 = reactionChoiceIndices2[r];
			if (choiceIndex1 != -1) {
				if (choiceIndex2 != -1) {
					if (choiceIndex1 == choiceIndex2)
						p *= (1 / 2.0) * x[choiceIndex1] * (x[choiceIndex1] - 1);
					else
						p *= x[choiceIndex1] * x[choiceIndex2];
				} else
					p *= x[choiceIndex1];
			}
			propensities[r] = p;
		}
	}

//	@Override
//	public void computePropensities(double t, double[] x, double[] propensities) {
//		// We don't check the length of x and propensities for performance reasons
//		for (int i = 0; i < stochasticRateParameters.length; i++) {
//			double p = stochasticRateParameters[i];
//			int[] choiceIndices = stochasticReactionChoiceIndices.get(i);
//			if (choiceIndices.length == 1)
//				p *= x[choiceIndices[0]];
//			else if (choiceIndices.length == 2)
//				if (choiceIndices[0] == choiceIndices[1])
//					p *= (1 / 2.0) * x[choiceIndices[0]] * (x[choiceIndices[1]] - 1);
//				else
//					p *= x[choiceIndices[0]] * x[choiceIndices[1]];
//			propensities[i] = p;
//		}
//		for (int i = stochasticRateParameters.length; i < stochasticRateParameters.length + hybridRateParameters.length; i++) {
//			double p = hybridRateParameters[i - stochasticRateParameters.length];
//			int[] choiceIndices = hybridReactionChoiceIndices.get(i - stochasticRateParameters.length);
//			if (choiceIndices.length == 1)
//				p *= x[choiceIndices[0]];
//			else if (choiceIndices.length == 2)
//				if (choiceIndices[0] == choiceIndices[1])
//					p *= (1 / 2.0) * x[choiceIndices[0]] * (x[choiceIndices[1]] - 1);
//				else
//					p *= x[choiceIndices[0]] * x[choiceIndices[1]];
//			propensities[i] = p;
//		}
//	}

	@Override
	public void updateState(int reaction, double t, double[] x) {
		// We don't check the length of x and the value of reaction for performance reasons
		double[] stochiometry = reactionStochiometries[reaction];
		for (int i = 0; i < stochiometry.length; i++) {
			x[i] += stochiometry[i];
		}
	}

//	@Override
//	public void updateState(int reaction, double t, double[] x) {
//		// We don't check the length of x and the value of reaction for performance reasons
//		int[] stochiometry;
//		if (reaction < stochasticRateParameters.length)
//			stochiometry = stochasticReactionStochiometries[reaction];
//		else
//			stochiometry = hybridReactionStochiometries[reaction];
//		for (int i = 0; i < stochiometry.length; i++) {
//			x[i] += stochiometry[i];
//		}
//	}

	@Override
	public int getNumberOfStates() {
		return getDimension();
	}

	@Override
	public boolean isTimeIndependent() {
		return true;
	}

}
