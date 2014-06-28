package ch.ethz.bhepp.hybridstochasticsimulation.math.distributions;

import org.apache.commons.math3.random.RandomDataGenerator;

public class PoissonDistribution implements UnivariateDistribution {

	public static long sample(RandomDataGenerator rdg, double lambda) {
		return rdg.nextPoisson(lambda);
	}

	private RandomDataGenerator rdg;
	private double lambda;

	public PoissonDistribution(RandomDataGenerator rdg, double lambda) {
		this.rdg = rdg;
		this.lambda = lambda;
	}

	@Override
	public double sample() {
		return PoissonDistribution.sample(rdg, lambda);
	}

	@Override
	public double getFirstMoment() {
		return lambda;
	}

	@Override
	public double getSecondMoment() {
		return lambda * lambda + lambda;
	}

}
