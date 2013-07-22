package ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders;

import javax.inject.Provider;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;

import com.google.inject.Inject;

public class MSHybridReactionNetworkProvider extends AbstractProvider<MSHybridReactionNetwork> {

	private Provider<UnaryBinaryReactionNetwork> networkProvider;

	@Inject
	public MSHybridReactionNetworkProvider(HierarchicalConfiguration config, Provider<UnaryBinaryReactionNetwork> networkProvider) {
		super(config, "ModelParameters");
		this.networkProvider = networkProvider;
	}

	@Override
	public MSHybridReactionNetwork get() {
		UnaryBinaryReactionNetwork network = networkProvider.get();
		double N = config().getDouble("N");
		double gamma = config().getDouble("gamma");
		double[] alpha = dataConfig().getDoubleArray("alpha", null);
//		if (alpha == null) {
//			alpha = new double[network.getNumberOfSpecies()];
//			Arrays.fill(alpha, 0.0);
//		}
		double[] beta = dataConfig().getDoubleArray("beta", null);
//		if (beta == null) {
//			beta = new double[network.getNumberOfReactions()];
//			Arrays.fill(beta, 0.0);
//		}
		return MSHybridReactionNetwork.createFrom(network, N, gamma, alpha, beta);
	}

}
