package ch.ethz.khammash.hybridstochasticsimulation.batch.providers;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.controllers.SimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.controllers.StochasticSimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.factories.RandomDataGeneratorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;

import com.google.inject.Inject;

public class StochasticSimulationControllerProvider extends AbstractSimulationControllerProvider<StochasticReactionNetworkModel> {

	@Inject
	public StochasticSimulationControllerProvider(HierarchicalConfiguration config, RandomDataGeneratorFactory rdgFactory) {
		super(config, rdgFactory);
	}

	@Override
	protected SimulationController<StochasticReactionNetworkModel> getController(int numOfThreads) {
		return new StochasticSimulationController(numOfThreads);
	}

}
