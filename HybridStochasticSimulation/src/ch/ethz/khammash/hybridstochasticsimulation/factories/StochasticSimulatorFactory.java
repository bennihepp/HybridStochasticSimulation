package ch.ethz.khammash.hybridstochasticsimulation.factories;

import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.StochasticSimulator;

public class StochasticSimulatorFactory
		implements SimulatorFactory<Simulator<StochasticReactionNetworkModel>> {

	@Override
	public Simulator<StochasticReactionNetworkModel> createSimulator(RandomDataGenerator rdg) {
		return new StochasticSimulator(rdg);
	}

}