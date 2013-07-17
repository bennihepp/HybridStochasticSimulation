package ch.ethz.khammash.hybridstochasticsimulation.math;

public interface MultivariateFunction {

	int getDimension();

	void computeValue(double[] x, double[] y);

}
