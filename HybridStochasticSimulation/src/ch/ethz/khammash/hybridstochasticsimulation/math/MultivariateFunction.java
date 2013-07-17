package ch.ethz.khammash.hybridstochasticsimulation.sandbox;

public interface MultivariateFunction {

	int getDimension();

	void computeValue(double[] x, double[] y);

}
