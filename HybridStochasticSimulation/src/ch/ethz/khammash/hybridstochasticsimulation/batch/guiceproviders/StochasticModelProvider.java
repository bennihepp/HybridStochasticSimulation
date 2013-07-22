package ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders;

import javax.inject.Provider;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.UnaryBinaryStochasticModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;

import com.google.inject.Inject;

public class StochasticModelProvider extends AbstractObjProvider<StochasticReactionNetworkModel> {

	private Provider<UnaryBinaryReactionNetwork> networkProvider;

	@Inject
	public StochasticModelProvider(HierarchicalConfiguration config, Provider<UnaryBinaryReactionNetwork> networkProvider) {
		super(config, "ModelParameters");
		this.networkProvider = networkProvider;
	}

	@Override
	public StochasticReactionNetworkModel get() {
		UnaryBinaryReactionNetwork network = networkProvider.get();
		return new UnaryBinaryStochasticModel(network);
	}

}
