package ch.ethz.khammash.hybridstochasticsimulation.batch.providers;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.factories.ModelFactory;
import ch.ethz.khammash.hybridstochasticsimulation.models.HybridModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModelAdapter;
import ch.ethz.khammash.hybridstochasticsimulation.models.UnaryBinaryDeterministicModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;

import com.google.inject.Inject;

public class DeterministicModelFactoryProvider extends AbstractProvider<ModelFactory<PDMPModel>> {

	private UnaryBinaryReactionNetwork network;

	@Inject
	public DeterministicModelFactoryProvider(HierarchicalConfiguration config, UnaryBinaryReactionNetwork network) {
		super(config, "ModelParameters");
		this.network = network;
	}

	@Override
	public ModelFactory<PDMPModel> get() {
		return new ModelFactory<PDMPModel>() {

			@Override
			public PDMPModel createModel() {
				UnaryBinaryDeterministicModel model = new UnaryBinaryDeterministicModel(network);
				return new PDMPModelAdapter<HybridModel>(model);
			}

		};
	}

}
