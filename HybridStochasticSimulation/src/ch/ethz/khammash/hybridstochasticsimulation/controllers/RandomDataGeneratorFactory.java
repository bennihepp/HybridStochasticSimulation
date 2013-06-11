package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import org.apache.commons.math3.random.RandomDataGenerator;

public interface RandomDataGeneratorFactory {
	public RandomDataGenerator createRandomDataGenerator();
}