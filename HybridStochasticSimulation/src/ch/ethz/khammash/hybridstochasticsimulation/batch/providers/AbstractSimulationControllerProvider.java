package ch.ethz.khammash.hybridstochasticsimulation.batch.providers;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.controllers.SimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.factories.RandomDataGeneratorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;

import com.google.inject.Inject;

public abstract class AbstractSimulationControllerProvider<T extends ReactionNetworkModel> extends AbstractProvider<SimulationController<T>> {

	private RandomDataGeneratorFactory rdgFactory;

	@Inject
	public AbstractSimulationControllerProvider(HierarchicalConfiguration config, RandomDataGeneratorFactory rdgFactory) {
		super(config, "SimulationParameters");
		this.rdgFactory = rdgFactory;
	}

	@Override
	public SimulationController<T> get() {
		int numOfThreads = config().getInt("numOfThreads");
		SimulationController<T> simCtrl = getController(numOfThreads);
		simCtrl.setRandomDataGeneratorFactory(rdgFactory);
		return simCtrl;
	}

	protected abstract SimulationController<T> getController(int numOfThreads);

}
