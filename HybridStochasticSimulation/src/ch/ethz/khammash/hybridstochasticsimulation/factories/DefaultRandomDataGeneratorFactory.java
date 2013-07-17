package ch.ethz.khammash.hybridstochasticsimulation.factories;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;

public class DefaultRandomDataGeneratorFactory implements RandomDataGeneratorFactory {

	private RandomGenerator baseRandomGenerator;
	private RandomDataGenerator baseRandomDataGenerator;

	public DefaultRandomDataGeneratorFactory() {
		baseRandomGenerator = new MersenneTwister();
	}

	public DefaultRandomDataGeneratorFactory(RandomGenerator rnd) {
		this.baseRandomGenerator = rnd;
	}

	public DefaultRandomDataGeneratorFactory(RandomDataGenerator rdg) {
		baseRandomDataGenerator = rdg;
	}

	@Override
	public RandomDataGenerator createRandomDataGenerator() {
		RandomDataGenerator rdg = new RandomDataGenerator();
		if (baseRandomGenerator != null) {
			long seed = baseRandomGenerator.nextLong();
			rdg.reSeed(seed);
		} else if (baseRandomDataGenerator != null) {
			long seed = baseRandomDataGenerator.nextLong(Long.MIN_VALUE, Long.MAX_VALUE);
			rdg.reSeed(seed);
		}
		return rdg;
	}

}