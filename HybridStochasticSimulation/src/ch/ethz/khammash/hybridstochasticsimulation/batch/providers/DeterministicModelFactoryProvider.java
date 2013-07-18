package ch.ethz.khammash.hybridstochasticsimulation.batch.providers;

import javax.inject.Provider;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.factories.ModelFactory;
import ch.ethz.khammash.hybridstochasticsimulation.models.HybridModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModelAdapter;
import ch.ethz.khammash.hybridstochasticsimulation.models.UnaryBinaryDeterministicModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;

import com.google.inject.Inject;

public class DeterministicModelFactoryProvider extends AbstractProvider<ModelFactory<PDMPModel>> {

	private Provider<UnaryBinaryReactionNetwork> networkProvider;

	@Inject
	public DeterministicModelFactoryProvider(HierarchicalConfiguration config, Provider<UnaryBinaryReactionNetwork> networkProvider) {
		super(config, "ModelParameters");
		this.networkProvider = networkProvider;
	}

	@Override
	public ModelFactory<PDMPModel> get() {
		return new ModelFactory<PDMPModel>() {

			@Override
			public PDMPModel createModel() {
				UnaryBinaryReactionNetwork network = networkProvider.get();
				UnaryBinaryDeterministicModel model = new UnaryBinaryDeterministicModel(network);
				return new PDMPModelAdapter<HybridModel>(model);
			}

		};
	}

}
