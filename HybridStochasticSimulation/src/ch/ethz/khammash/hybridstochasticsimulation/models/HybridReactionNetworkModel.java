package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.Arrays;
import java.util.LinkedList;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import ch.ethz.khammash.hybridstochasticsimulation.networks.HybridReactionNetwork;


public class HybridReactionNetworkModel implements HybridModel,
		FirstOrderDifferentialEquations, StochasticReactionNetworkModel {

	private int dimension;

	private int[] stochasticReactionIndices;
	private int[] deterministicReactionIndices;

	private double[] modelRateParameters;
	private int[] reactionChoiceIndices1;
	private int[] reactionChoiceIndices2;
	private double[][] reactionStochiometries;

	public HybridReactionNetworkModel(HybridReactionNetwork hrn) {
		dimension = hrn.getNumberOfSpecies();
		modelRateParameters = new double[hrn.getNumberOfReactions()];
		reactionChoiceIndices1 = new int[hrn.getNumberOfReactions()];
		reactionChoiceIndices2 = new int[hrn.getNumberOfReactions()];
		reactionStochiometries = new double[hrn.getNumberOfReactions()][hrn.getNumberOfSpecies()];

		init(hrn);
	}

	final private void init(HybridReactionNetwork hrn) {

		LinkedList<Integer> stochasticReactionIndicesList = new LinkedList<Integer>();
		LinkedList<Integer> deterministicReactionIndicesList = new LinkedList<Integer>();

		for (int r = 0; r < hrn.getNumberOfReactions(); r++) {
			for (int s = 0; s < hrn.getNumberOfSpecies(); s++) {
				reactionStochiometries[r][s] = hrn.getStochiometry(s, r);
			}
			modelRateParameters[r] = hrn.getRateParameter(r);
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
				if (choiceIndices[0] == choiceIndices[1])
					modelRateParameters[r] *= 1 / 2.0;
				break;
			}
			if (choiceIndices.length == 2 && choiceIndices[0] == choiceIndices[1])
				modelRateParameters[r] /= 2.0;
			HybridReactionNetwork.ReactionType rt = hrn.getReactionType(r);
			switch (rt) {
			case STOCHASTIC:
				stochasticReactionIndicesList.add(Integer.valueOf(r));
				break;
			case DETERMINISTIC:
				deterministicReactionIndicesList.add(Integer.valueOf(r));
				break;
			}
		}

		stochasticReactionIndices = ArrayUtils.toPrimitive(stochasticReactionIndicesList.toArray(new Integer[0]));
		deterministicReactionIndices = ArrayUtils.toPrimitive(deterministicReactionIndicesList.toArray(new Integer[0]));
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
			double v = modelRateParameters[r];
			int choiceIndex1 = reactionChoiceIndices1[r];
			int choiceIndex2 = reactionChoiceIndices2[r];
			if (choiceIndex1 != -1) {
				v *= x[choiceIndex1];
				if (choiceIndex2 != -1) {
					double q = x[choiceIndex2];
					if (choiceIndex1 == choiceIndex2)
						q -= 1;
					v *= q;
				}
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
		return modelRateParameters.length;
	}

	@Override
	final public double computePropensity(int reaction, double t, double[] x) {
		double p = modelRateParameters[reaction];
		int choiceIndex1 = reactionChoiceIndices1[reaction];
		int choiceIndex2 = reactionChoiceIndices2[reaction];
		if (choiceIndex1 != -1) {
			p *= x[choiceIndex1];
			if (choiceIndex2 != -1) {
				double q = x[choiceIndex2];
				if (choiceIndex1 == choiceIndex2)
					q -= 1;
				p *= q;
			}
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
