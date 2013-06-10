package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.List;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

public class DeterministicModel implements FirstOrderDifferentialEquations {

	private int dimension;
	private double[] rateParameters;
	private int[] reactionChoiceIndices1;
	private int[] reactionChoiceIndices2;
	private double[][] reactionStochiometries;

	public DeterministicModel(ReactionNetwork net) {
		dimension = net.getNumberOfSpecies();
		rateParameters = new double[net.getNumberOfReactions()];
		reactionChoiceIndices1 = new int[net.getNumberOfReactions()];
		reactionChoiceIndices2 = new int[net.getNumberOfReactions()];
		reactionStochiometries = new double[net.getNumberOfReactions()][net.getNumberOfSpecies()];
		init(net);
	}

	final protected void init(ReactionNetwork net) {
		List<int[]> choiceIndicesList = net.getChoiceIndicesList();
		for (int r = 0; r < net.getNumberOfReactions(); r++) {
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
			rateParameters[r] = net.getRateParameter(r);
			for (int s = 0; s < net.getNumberOfSpecies(); s++)
				reactionStochiometries[r][s] = net.getStochiometry(s, r);
		}
	}

	public FirstOrderDifferentialEquations getFirstOrderDifferentialEquations() {
		return this;
	}

	@Override
	public int getDimension() {
		return dimension;
	}

	@Override
	public void computeDerivatives(double t, double[] x, double[] xDot)
			throws MaxCountExceededException, DimensionMismatchException {
		// We don't check the length of x and propensities for performance reasons
		for (int s = 0; s < xDot.length; s++)
			xDot[s] = 0;
		for (int r=0; r < rateParameters.length; r++) {
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

}
