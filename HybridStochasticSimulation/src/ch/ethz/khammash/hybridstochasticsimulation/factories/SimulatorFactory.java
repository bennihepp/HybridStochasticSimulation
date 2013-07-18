package ch.ethz.khammash.hybridstochasticsimulation.factories;

import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;

public interface SimulatorFactory<T extends ReactionNetworkModel> {

	Simulator<T> createSimulator(RandomDataGenerator rdg);

}
