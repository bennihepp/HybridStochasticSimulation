package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;

public interface SimulatorFactory<S extends Simulator<?,?>> {

	public S createSimulator(RandomDataGenerator rdg);

}
