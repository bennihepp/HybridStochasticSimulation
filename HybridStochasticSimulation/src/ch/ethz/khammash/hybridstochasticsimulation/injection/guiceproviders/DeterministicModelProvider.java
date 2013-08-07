package ch.ethz.khammash.hybridstochasticsimulation.injection.guiceproviders;

import javax.inject.Provider;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.models.HybridModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModelAdapter;
import ch.ethz.khammash.hybridstochasticsimulation.models.UnaryBinaryDeterministicModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;

import com.google.inject.Inject;

public class DeterministicModelProvider extends AbstractObjProvider<PDMPModel> {

	private Provider<UnaryBinaryReactionNetwork> networkProvider;

	@Inject
	public DeterministicModelProvider(HierarchicalConfiguration config, Provider<UnaryBinaryReactionNetwork> networkProvider) {
		super(config, "ModelParameters");
		this.networkProvider = networkProvider;
	}

	@Override
	public PDMPModel get() {
		UnaryBinaryReactionNetwork network = networkProvider.get();
		UnaryBinaryDeterministicModel model = new UnaryBinaryDeterministicModel(network);
		return new PDMPModelAdapter<HybridModel>(model);
	}

}
