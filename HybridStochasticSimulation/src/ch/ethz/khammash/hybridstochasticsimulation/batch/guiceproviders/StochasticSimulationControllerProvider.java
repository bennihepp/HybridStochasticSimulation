package ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders;

import org.apache.commons.configuration.HierarchicalConfiguration;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.controllers.SimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.controllers.StochasticSimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.providers.ObjProvider;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;

import com.google.inject.Inject;

public class StochasticSimulationControllerProvider extends AbstractSimulationControllerProvider<StochasticReactionNetworkModel> {

	private ObjProvider<Simulator<StochasticReactionNetworkModel>> simulatorProvider;

	@Inject
	public StochasticSimulationControllerProvider(HierarchicalConfiguration config,
			ObjProvider<Simulator<StochasticReactionNetworkModel>> simulatorProvider,
			ObjProvider<RandomDataGenerator> rdgProvider) {
		super(config, rdgProvider);
		this.simulatorProvider = simulatorProvider;
	}

	@Override
	protected SimulationController<StochasticReactionNetworkModel> getController(int numOfThreads) {
		StochasticSimulationController simCtrl = new StochasticSimulationController(numOfThreads);
		simCtrl.setSimulatorProvider(simulatorProvider);
		return simCtrl;
	}

}
