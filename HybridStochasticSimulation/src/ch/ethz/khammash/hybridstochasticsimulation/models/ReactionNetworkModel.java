package ch.ethz.khammash.hybridstochasticsimulation.models;

public interface ReactionNetworkModel {

	public int getStateDimension();

	public int getPropensityDimension();

	public void computePropensities(double t, double[] x, double[] propensities);

	public void updateState(int reaction, double t, double[] x);

}
