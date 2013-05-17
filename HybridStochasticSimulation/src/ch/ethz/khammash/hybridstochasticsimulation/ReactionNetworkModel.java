package ch.ethz.khammash.hybridstochasticsimulation;

public interface ReactionNetworkModel {

	public int getNumberOfSpecies();

	public int getPropensityDimension();

	public void computePropensities(double t, double[] x, double[] propensities);

	public void updateState(int reaction, double t, double[] x);

}
