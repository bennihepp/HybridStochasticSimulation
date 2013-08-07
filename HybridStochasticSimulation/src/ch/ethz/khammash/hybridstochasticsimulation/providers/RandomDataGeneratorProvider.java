package ch.ethz.khammash.hybridstochasticsimulation.providers;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;

public class RandomDataGeneratorProvider implements ObjProvider<RandomDataGenerator> {

	private RandomGenerator baseRandomGenerator;
	private RandomDataGenerator baseRandomDataGenerator;

	public RandomDataGeneratorProvider() {
		baseRandomGenerator = new MersenneTwister();
	}

	public RandomDataGeneratorProvider(long seed) {
		baseRandomGenerator = new MersenneTwister(seed);
	}

	public RandomDataGeneratorProvider(RandomGenerator rnd) {
		this.baseRandomGenerator = rnd;
	}

	public RandomDataGeneratorProvider(RandomDataGenerator rdg) {
		baseRandomDataGenerator = rdg;
	}

	@Override
	public RandomDataGenerator get() {
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