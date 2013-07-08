package ch.ethz.khammash.hybridstochasticsimulation.factories;

import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;

public interface SimulatorFactory<S extends Simulator<?>> {

	S createSimulator(RandomDataGenerator rdg);

}
