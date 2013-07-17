package ch.ethz.khammash.hybridstochasticsimulation.batch.providers;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;

import com.google.inject.Inject;

public class MSHybridReactionNetworkProvider extends AbstractProvider<MSHybridReactionNetwork> {

	private UnaryBinaryReactionNetwork network;

	@Inject
	public MSHybridReactionNetworkProvider(HierarchicalConfiguration config, UnaryBinaryReactionNetwork network) {
		super(config, "ModelParameters");
		this.network = network;
	}

	@Override
	public MSHybridReactionNetwork get() {
		double N = config().getDouble("N");
		double gamma = config().getDouble("gamma");
		double[] alpha = dataConfig().getDoubleArray("alpha");
		double[] beta = dataConfig().getDoubleArray("beta");
		return MSHybridReactionNetwork.createFrom(network, N, gamma, alpha, beta);
	}

}
