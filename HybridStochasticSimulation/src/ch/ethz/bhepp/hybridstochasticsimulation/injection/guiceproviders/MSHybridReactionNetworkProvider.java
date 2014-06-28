package ch.ethz.bhepp.hybridstochasticsimulation.injection.guiceproviders;

import javax.inject.Provider;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.bhepp.hybridstochasticsimulation.networks.MSHybridReactionNetwork;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.MassActionReactionNetwork;

import com.google.inject.Inject;

public class MSHybridReactionNetworkProvider extends AbstractProvider<MSHybridReactionNetwork> {

	private Provider<MassActionReactionNetwork> networkProvider;

	@Inject
	public MSHybridReactionNetworkProvider(HierarchicalConfiguration config, Provider<MassActionReactionNetwork> networkProvider) {
		super(config, "ModelParameters");
		this.networkProvider = networkProvider;
	}

	@Override
	public MSHybridReactionNetwork get() {
		MassActionReactionNetwork network = networkProvider.get();
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
