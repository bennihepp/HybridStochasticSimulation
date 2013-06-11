package ch.ethz.khammash.hybridstochasticsimulation.models;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;


public interface HybridModel {

	public FirstOrderDifferentialEquations getDeterministicModel();

	public StochasticReactionNetworkModel getStochasticModel();

	public boolean hasDeterministicPart();

	public boolean isTimeIndependent();

}
