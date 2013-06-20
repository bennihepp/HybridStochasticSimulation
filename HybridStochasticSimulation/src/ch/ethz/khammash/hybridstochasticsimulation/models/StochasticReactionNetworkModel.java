package ch.ethz.khammash.hybridstochasticsimulation.models;

public interface StochasticReactionNetworkModel extends ReactionNetworkModel {

	double computePropensity(int reaction, double t, double[] x);

	void computePropensities(double t, double[] x, double[] propensities);

	void updateState(int reaction, double t, double[] x);

	boolean isTimeIndependent();

}
