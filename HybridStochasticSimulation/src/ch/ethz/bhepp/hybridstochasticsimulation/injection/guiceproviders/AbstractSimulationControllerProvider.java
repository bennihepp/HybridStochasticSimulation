package ch.ethz.bhepp.hybridstochasticsimulation.injection.guiceproviders;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.bhepp.hybridstochasticsimulation.controllers.SimulationController;
import ch.ethz.bhepp.hybridstochasticsimulation.models.ReactionNetworkModel;

import com.google.inject.Inject;

public abstract class AbstractSimulationControllerProvider<T extends ReactionNetworkModel> extends AbstractProvider<SimulationController<T>> {

	@Inject
	public AbstractSimulationControllerProvider(HierarchicalConfiguration config) {
		super(config, "SimulationParameters");
	}

	@Override
	public SimulationController<T> get() {
		int numOfThreads = config().getInt("numOfThreads");
		SimulationController<T> simCtrl = getController(numOfThreads);
		return simCtrl;
	}

	protected abstract SimulationController<T> getController(int numOfThreads);

}
