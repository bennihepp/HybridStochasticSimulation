package ch.ethz.khammash.hybridstochasticsimulation.factories;

import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;

public interface SimulatorFactory<S extends Simulator<?,?>> {

	// TODO: Create without RandomDataGenerator
	S createSimulator(RandomDataGenerator rdg);

}
