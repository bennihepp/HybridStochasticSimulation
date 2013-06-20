package ch.ethz.khammash.hybridstochasticsimulation.factories;

import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.StochasticSimulator;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;

public class StochasticSimulatorFactory<T extends StochasticReactionNetworkModel>
		implements SimulatorFactory<Simulator<T, TrajectoryRecorder<T>>> {

	@Override
	public StochasticSimulator<T> createSimulator(RandomDataGenerator rdg) {
		return new StochasticSimulator<T>(rdg);
	}

}