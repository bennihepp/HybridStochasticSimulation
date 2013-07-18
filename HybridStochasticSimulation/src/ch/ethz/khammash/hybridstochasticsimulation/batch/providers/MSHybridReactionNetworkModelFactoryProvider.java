package ch.ethz.khammash.hybridstochasticsimulation.batch.providers;

import javax.inject.Provider;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.factories.ModelFactory;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPMSHRNModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork;

import com.google.inject.Inject;

public class MSHybridReactionNetworkModelFactoryProvider extends AbstractProvider<ModelFactory<PDMPModel>> {

	private Provider<MSHybridReactionNetwork> hrnProvider;

	@Inject
	public MSHybridReactionNetworkModelFactoryProvider(HierarchicalConfiguration config, Provider<MSHybridReactionNetwork> hrnProvider) {
		super(config, "ModelParameters");
		this.hrnProvider = hrnProvider;
	}

	@Override
	public ModelFactory<PDMPModel> get() {
		return new ModelFactory<PDMPModel>() {

			@Override
			public PDMPModel createModel() {
				MSHybridReactionNetwork hrn = hrnProvider.get();
				PDMPMSHRNModel hrnModel = new PDMPMSHRNModel(hrn);
				hrnModel.setExposeOptionalState(config().getBoolean("exposeOptionalState", false));
				return hrnModel;
			}

		};
	}

}
