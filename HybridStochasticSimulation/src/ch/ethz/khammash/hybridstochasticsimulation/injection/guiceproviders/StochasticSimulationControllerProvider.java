package ch.ethz.khammash.hybridstochasticsimulation.injection.guiceproviders;

import org.apache.commons.configuration.HierarchicalConfiguration;

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
			ObjProvider<Simulator<StochasticReactionNetworkModel>> simulatorProvider) {
		super(config);
		this.simulatorProvider = simulatorProvider;
	}

	@Override
	protected SimulationController<StochasticReactionNetworkModel> getController(int numOfThreads) {
		StochasticSimulationController simCtrl = new StochasticSimulationController(simulatorProvider, numOfThreads);
		return simCtrl;
	}

}
