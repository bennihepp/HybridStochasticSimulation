package ch.ethz.bhepp.hybridstochasticsimulation.math;

public interface MultivariateFunction {

	int getDimension();

	void computeValue(double[] x, double[] y);

}
