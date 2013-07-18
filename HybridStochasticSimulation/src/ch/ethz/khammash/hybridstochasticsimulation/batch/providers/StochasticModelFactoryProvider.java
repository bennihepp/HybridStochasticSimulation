package ch.ethz.khammash.hybridstochasticsimulation.batch.providers;

import javax.inject.Provider;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.factories.ModelFactory;
import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.UnaryBinaryStochasticModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;

import com.google.inject.Inject;

public class StochasticModelFactoryProvider extends AbstractProvider<ModelFactory<StochasticReactionNetworkModel>> {

	private Provider<UnaryBinaryReactionNetwork> networkProvider;

	@Inject
	public StochasticModelFactoryProvider(HierarchicalConfiguration config, Provider<UnaryBinaryReactionNetwork> networkProvider) {
		super(config, "ModelParameters");
		this.networkProvider = networkProvider;
	}

	@Override
	public ModelFactory<StochasticReactionNetworkModel> get() {
		return new ModelFactory<StochasticReactionNetworkModel>() {

			@Override
			public StochasticReactionNetworkModel createModel() {
				UnaryBinaryReactionNetwork network = networkProvider.get();
				return new UnaryBinaryStochasticModel(network);
			}

		};
	}

}
