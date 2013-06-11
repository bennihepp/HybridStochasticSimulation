package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.Arrays;
import java.util.LinkedList;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import ch.ethz.khammash.hybridstochasticsimulation.networks.HybridReactionNetwork;


// TODO: Group computation of choices for each reaction

public class HybridReactionNetworkModel implements HybridModel,
		FirstOrderDifferentialEquations, StochasticReactionNetworkModel {

	private int dimension;

	private int[] stochasticReactionIndices;
	private int[] deterministicReactionIndices;

	private double[] rateParameters;
	private int[] reactionChoiceIndices1;
	private int[] reactionChoiceIndices2;
	private double[][] reactionStochiometries;

	public HybridReactionNetworkModel(HybridReactionNetwork hrn) {
		dimension = hrn.getNumberOfSpecies();
		rateParameters = new double[hrn.getNumberOfReactions()];
		reactionChoiceIndices1 = new int[hrn.getNumberOfReactions()];
		reactionChoiceIndices2 = new int[hrn.getNumberOfReactions()];
		reactionStochiometries = new double[hrn.getNumberOfReactions()][hrn.getNumberOfSpecies()];

		init(hrn);
	}

	final private void init(HybridReactionNetwork hrn) {

		LinkedList<Integer> stochasticReactionIndicesList = new LinkedList<Integer>();
		LinkedList<Integer> deterministicReactionIndicesList = new LinkedList<Integer>();

		for (int r = 0; r < hrn.getNumberOfReactions(); r++) {
			int[] choiceIndices = hrn.getChoiceIndices(r);
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

	@Override
	public FirstOrderDifferentialEquations getDeterministicModel() {
		return this;
	}

	@Override
	public StochasticReactionNetworkModel getStochasticModel() {
		return this;
	}

	@Override
	public boolean hasDeterministicPart() {
		return deterministicReactionIndices.length > 0;
	}

	@Override
	public boolean isTimeIndependent() {
		return true;
	}

	@Override
	public int getDimension() {
		return dimension;
	}

	@Override
	public void computeDerivatives(double t, double[] x, double[] xDot) {
		// We don't check the length of x and propensities for performance reasons
		Arrays.fill(xDot, 0, getNumberOfSpecies(), 0.0);
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

	@Override
	public int getNumberOfSpecies() {
		return getDimension();
	}

	@Override
	public int getNumberOfReactions() {
		return rateParameters.length;
	}

	@Override
	final public double computePropensity(int reaction, double t, double[] x) {
		double p = rateParameters[reaction];
		int choiceIndex1 = reactionChoiceIndices1[reaction];
		int choiceIndex2 = reactionChoiceIndices2[reaction];
		if (choiceIndex1 != -1) {
			if (choiceIndex2 != -1) {
				if (choiceIndex1 == choiceIndex2)
					p *= (1 / 2.0) * x[choiceIndex1] * (x[choiceIndex1] - 1);
				else
					p *= x[choiceIndex1] * x[choiceIndex2];
			} else
				p *= x[choiceIndex1];
		}
		return p;
	}

	@Override
	public void computePropensities(double t, double[] x, double[] propensities) {
		// We don't check the length of x and propensities for performance reasons
		Arrays.fill(propensities, 0, getNumberOfReactions(), 0.0);
		for (int r : stochasticReactionIndices)
			propensities[r] = computePropensity(r, t, x);
	}

	@Override
	public void updateState(int reaction, double t, double[] x) {
		// We don't check the length of x and the value of reaction for performance reasons
		double[] stochiometry = reactionStochiometries[reaction];
		for (int i = 0; i < stochiometry.length; i++) {
			x[i] += stochiometry[i];
		}
	}

}
