package ch.ethz.khammash.hybridstochasticsimulation.batch.providers;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.factories.ModelFactory;
import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.UnaryBinaryStochasticModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;

import com.google.inject.Inject;

public class StochasticModelFactoryProvider extends AbstractProvider<ModelFactory<StochasticReactionNetworkModel>> {

	private UnaryBinaryReactionNetwork network;

	@Inject
	public StochasticModelFactoryProvider(HierarchicalConfiguration config, UnaryBinaryReactionNetwork network) {
		super(config, "ModelParameters");
		this.network = network;
	}

	@Override
	public ModelFactory<StochasticReactionNetworkModel> get() {
		return new ModelFactory<StochasticReactionNetworkModel>() {

			@Override
			public StochasticReactionNetworkModel createModel() {
				return new UnaryBinaryStochasticModel(network);
			}

		};
	}

}
