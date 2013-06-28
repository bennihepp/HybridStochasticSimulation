package ch.ethz.khammash.hybridstochasticsimulation.models;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;


public interface HybridModel extends ReactionNetworkModel {

	FirstOrderDifferentialEquations getVectorField();

	StochasticReactionNetworkModel getTransitionMeasure();

	void computeDerivativesAndPropensities(double t, double[] x, double[] xDot, double[] propensities);

	double computeDerivativesAndPropensitiesSum(double t, double[] x, double[] xDot);

	boolean hasDeterministicPart();

	boolean isTimeIndependent();

}
