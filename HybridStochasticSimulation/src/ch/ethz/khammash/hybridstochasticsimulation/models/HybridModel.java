package ch.ethz.khammash.hybridstochasticsimulation.models;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;


public interface HybridModel extends FirstOrderDifferentialEquations, ReactionNetworkModel {

	public boolean hasDeterministicPart();

	public boolean isTimeIndependent();

}
