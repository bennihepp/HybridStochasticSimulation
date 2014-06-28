package ch.ethz.bhepp.hybridstochasticsimulation.math.distributions;

public interface MultivariateDistribution {

	double[] sample();

	void sample(double[] y);

	double[] getFirstMoment();

	void getFirstMoment(double[] firstMoment);

	double[][] getSecondMoment();

	void getSecondMoment(double[][] secondMoment);

}
