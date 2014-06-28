package ch.ethz.bhepp.hybridstochasticsimulation.injection.guiceproviders;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.bhepp.hybridstochasticsimulation.controllers.SimulationController;
import ch.ethz.bhepp.hybridstochasticsimulation.controllers.StochasticSimulationController;
import ch.ethz.bhepp.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.bhepp.hybridstochasticsimulation.providers.ObjProvider;
import ch.ethz.bhepp.hybridstochasticsimulation.simulators.Simulator;

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
