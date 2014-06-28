package ch.ethz.bhepp.hybridstochasticsimulation.math.distributions;

public interface UnivariateDistribution {

	double sample();

	double getFirstMoment();

	double getSecondMoment();

}
