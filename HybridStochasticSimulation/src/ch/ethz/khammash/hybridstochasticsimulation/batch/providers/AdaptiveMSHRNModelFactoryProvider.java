package ch.ethz.khammash.hybridstochasticsimulation.batch.providers;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.batch.providers.AbstractProvider;
import ch.ethz.khammash.hybridstochasticsimulation.factories.ModelFactory;
import ch.ethz.khammash.hybridstochasticsimulation.models.AdaptiveMSHRNModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.AdaptiveMSHRN;

import com.google.inject.Inject;

public class AdaptiveMSHRNModelFactoryProvider extends AbstractProvider<ModelFactory<PDMPModel>> {

	private AdaptiveMSHRN hrn;

	@Inject
	public AdaptiveMSHRNModelFactoryProvider(HierarchicalConfiguration config, AdaptiveMSHRN hrn) {
		super(config, "ModelParameters");
		this.hrn = hrn;
	}

	@Override
	public ModelFactory<PDMPModel> get() {
		return new ModelFactory<PDMPModel>() {

			@Override
			public PDMPModel createModel() {
				AdaptiveMSHRNModel hrnModel = new AdaptiveMSHRNModel(hrn);
				hrnModel.setExposeOptionalState(config().getBoolean("exposeOptionalState", false));
				return hrnModel;
			}

		};
	}

}
