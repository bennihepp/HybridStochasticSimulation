package ch.ethz.khammash.hybridstochasticsimulation.models;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import ch.ethz.khammash.hybridstochasticsimulation.networks.ReactionNetwork;


public class HybridModelAdapter implements HybridModel {

	private FirstOrderDifferentialEquations vectorField;
	private StochasticReactionNetworkModel transitionMeasure;
	private double[] propTmpVector;

	public HybridModelAdapter(FirstOrderDifferentialEquations deterministicModel, StochasticReactionNetworkModel stochasticModel) {
		if (stochasticModel.getNumberOfSpecies() != deterministicModel.getDimension())
			throw new UnsupportedOperationException("Expected stochasticModel.getNumberOfSpecies() == deterministicModel.getDimension()");
		this.vectorField = deterministicModel;
		this.transitionMeasure = stochasticModel;
		propTmpVector = new double[stochasticModel.getNumberOfReactions()];
	}

	@Override
	public int getNumberOfSpecies() {
		return transitionMeasure.getNumberOfSpecies();
	}

	@Override
	public int getNumberOfReactions() {
		return transitionMeasure.getNumberOfReactions();
	}

	@Override
	public boolean hasDeterministicPart() {
		return true;
	}

	@Override
	public boolean isTimeIndependent() {
		return false;
	}

	@Override
	public FirstOrderDifferentialEquations getVectorField() {
		return vectorField;
	}

	@Override
	public StochasticReactionNetworkModel getTransitionMeasure() {
		return transitionMeasure;
	}

	@Override
	public void computeDerivativesAndPropensities(double t, double[] x,
			double[] xDot, double[] propensities) {
		getTransitionMeasure().computePropensities(t, x, propensities);
		getVectorField().computeDerivatives(t, x, xDot);
	}

	@Override
	public double computeDerivativesAndPropensitiesSum(double t, double[] x,  double[] xDot) {
		computeDerivativesAndPropensities(t, x, xDot, propTmpVector);
		double propSum = 0.0;
		for (int r=0; r < propTmpVector.length; r++)
			propSum += propTmpVector[r];
		return propSum;
	}

	@Override
	public ReactionNetwork getNetwork() {
		return transitionMeasure.getNetwork();
	}

}
