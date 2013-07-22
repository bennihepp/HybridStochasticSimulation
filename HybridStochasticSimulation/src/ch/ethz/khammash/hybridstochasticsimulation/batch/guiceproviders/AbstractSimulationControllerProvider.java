package ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders;

import org.apache.commons.configuration.HierarchicalConfiguration;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.controllers.SimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.providers.ObjProvider;

import com.google.inject.Inject;

public abstract class AbstractSimulationControllerProvider<T extends ReactionNetworkModel> extends AbstractProvider<SimulationController<T>> {

	private ObjProvider<RandomDataGenerator> rdgProvider;

	@Inject
	public AbstractSimulationControllerProvider(HierarchicalConfiguration config, ObjProvider<RandomDataGenerator> rdgProvider) {
		super(config, "SimulationParameters");
		this.rdgProvider = rdgProvider;
	}

	@Override
	public SimulationController<T> get() {
		int numOfThreads = config().getInt("numOfThreads");
		SimulationController<T> simCtrl = getController(numOfThreads);
		simCtrl.setRandomDataGeneratorProvider(rdgProvider);
		return simCtrl;
	}

	protected abstract SimulationController<T> getController(int numOfThreads);

}
