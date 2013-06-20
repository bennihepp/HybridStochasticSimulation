package ch.ethz.khammash.hybridstochasticsimulation.models;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;


public interface HybridModel extends ReactionNetworkModel {

	FirstOrderDifferentialEquations getDeterministicModel();

	StochasticReactionNetworkModel getStochasticModel();

	boolean hasDeterministicPart();

	boolean isTimeIndependent();

}
