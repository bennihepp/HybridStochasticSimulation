package ch.ethz.bhepp.hybridstochasticsimulation.injection.guiceproviders;

import javax.inject.Provider;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.bhepp.hybridstochasticsimulation.models.HybridModel;
import ch.ethz.bhepp.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.bhepp.hybridstochasticsimulation.models.PDMPModelAdapter;
import ch.ethz.bhepp.hybridstochasticsimulation.models.UnaryBinaryDeterministicModel;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.MassActionReactionNetwork;

import com.google.inject.Inject;

public class DeterministicModelProvider extends AbstractObjProvider<PDMPModel> {

	private Provider<MassActionReactionNetwork> networkProvider;

	@Inject
	public DeterministicModelProvider(HierarchicalConfiguration config, Provider<MassActionReactionNetwork> networkProvider) {
		super(config, "ModelParameters");
		this.networkProvider = networkProvider;
	}

	@Override
	public PDMPModel get() {
		MassActionReactionNetwork network = networkProvider.get();
		UnaryBinaryDeterministicModel model = new UnaryBinaryDeterministicModel(network);
		return new PDMPModelAdapter<HybridModel>(model);
	}

}
