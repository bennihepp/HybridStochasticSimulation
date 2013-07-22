package ch.ethz.khammash.hybridstochasticsimulation.providers;

import org.apache.commons.math3.random.RandomDataGenerator;

import com.google.inject.Inject;

import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.StochasticSimulator;

public class StochasticSimulatorProvider
		implements ObjProvider<Simulator<StochasticReactionNetworkModel>> {

	private ObjProvider<RandomDataGenerator> rdgProvider;

	@Inject
	public StochasticSimulatorProvider(ObjProvider<RandomDataGenerator> rdgProvider) {
		this.rdgProvider = rdgProvider;
	}

	@Override
	public Simulator<StochasticReactionNetworkModel> get() {
		return new StochasticSimulator(rdgProvider.get());
	}

}