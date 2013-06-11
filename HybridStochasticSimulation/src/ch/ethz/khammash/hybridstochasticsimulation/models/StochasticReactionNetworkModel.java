package ch.ethz.khammash.hybridstochasticsimulation.models;

public interface StochasticReactionNetworkModel {

	public int getNumberOfSpecies();

	public int getNumberOfReactions();

	public double computePropensity(int reaction, double t, double[] x);

	public void computePropensities(double t, double[] x, double[] propensities);

	public void updateState(int reaction, double t, double[] x);

	public boolean isTimeIndependent();

}
