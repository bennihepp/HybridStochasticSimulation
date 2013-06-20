package ch.ethz.khammash.hybridstochasticsimulation.factories;

import org.apache.commons.math3.random.RandomDataGenerator;

public interface RandomDataGeneratorFactory {

	RandomDataGenerator createRandomDataGenerator();

}
