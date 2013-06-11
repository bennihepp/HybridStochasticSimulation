package ch.ethz.khammash.hybridstochasticsimulation.models;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

public interface DeterministicReactionNetworkModel extends FirstOrderDifferentialEquations {

	public int getNumberOfSpecies();

	public int getNumberOfReactions();

	public boolean isTimeIndependent();

}
