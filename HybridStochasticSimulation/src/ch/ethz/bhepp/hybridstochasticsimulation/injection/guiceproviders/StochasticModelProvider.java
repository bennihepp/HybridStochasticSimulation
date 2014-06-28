package ch.ethz.bhepp.hybridstochasticsimulation.injection.guiceproviders;

import javax.inject.Provider;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.bhepp.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.bhepp.hybridstochasticsimulation.models.UnaryBinaryStochasticModel;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.MassActionReactionNetwork;

import com.google.inject.Inject;

public class StochasticModelProvider extends AbstractObjProvider<StochasticReactionNetworkModel> {

	private Provider<MassActionReactionNetwork> networkProvider;

	@Inject
	public StochasticModelProvider(HierarchicalConfiguration config, Provider<MassActionReactionNetwork> networkProvider) {
		super(config, "ModelParameters");
		this.networkProvider = networkProvider;
	}

	@Override
	public StochasticReactionNetworkModel get() {
		MassActionReactionNetwork network = networkProvider.get();
		return new UnaryBinaryStochasticModel(network);
	}

}
