package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import ch.ethz.khammash.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;

public class UnaryBinaryDeterministicModel implements HybridModel, FirstOrderDifferentialEquations, StochasticReactionNetworkModel {

	private int dimension;
	private double[] rateParameters;
	private int[] reactionChoiceIndices1;
	private int[] reactionChoiceIndices2;
	private double[][] reactionStochiometries;
	private double[] propTmpVector;

	public UnaryBinaryDeterministicModel(DefaultUnaryBinaryReactionNetwork net) {
		dimension = net.getNumberOfSpecies();
		rateParameters = new double[net.getNumberOfReactions()];
		reactionChoiceIndices1 = new int[net.getNumberOfReactions()];
		reactionChoiceIndices2 = new int[net.getNumberOfReactions()];
		reactionStochiometries = new double[net.getNumberOfReactions()][net.getNumberOfSpecies()];
		propTmpVector = new double[net.getNumberOfReactions()];
		init(net);
	}

	final private void init(DefaultUnaryBinaryReactionNetwork net) {
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

	@Override
	public int getNumberOfSpecies() {
		return dimension;
	}

	@Override
	public int getNumberOfReactions() {
		return rateParameters.length;
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

	@Override
	public boolean isTimeIndependent() {
		return true;
	}

	@Override
	public FirstOrderDifferentialEquations getVectorField() {
		return this;
	}

	@Override
	public StochasticReactionNetworkModel getTransitionMeasure() {
		return this;
	}

	@Override
	public boolean hasDeterministicPart() {
		return true;
	}

	@Override
	public double computePropensity(int reaction, double t, double[] x) {
		return 0;
	}

	@Override
	public void computePropensities(double t, double[] x, double[] propensities) {
		Arrays.fill(propensities, 0, getNumberOfReactions(), 0.0);
	}

	@Override
	public void changeState(int reaction, double t, double[] x) {
	}

	@Override
	public void computeDerivativesAndPropensities(double t, double[] x, double[] xDot, double[] propensities) {
		computePropensities(t, x, propensities);
		computeDerivatives(t, x, xDot);
	}

	@Override
	public double computeDerivativesAndPropensitiesSum(double t, double[] x, double[] xDot) {
		computeDerivativesAndPropensities(t, x, xDot, propTmpVector);
		double propSum = 0.0;
		for (int r=0; r < propTmpVector.length; r++)
			propSum += propTmpVector[r];
		return propSum;
	}

}
