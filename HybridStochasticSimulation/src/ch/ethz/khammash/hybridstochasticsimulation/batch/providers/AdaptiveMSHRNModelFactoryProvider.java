package ch.ethz.khammash.hybridstochasticsimulation.batch.providers;

import javax.inject.Provider;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.factories.ModelFactory;
import ch.ethz.khammash.hybridstochasticsimulation.models.AdaptiveMSHRNModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.AdaptiveMSHRN;

import com.google.inject.Inject;

public class AdaptiveMSHRNModelFactoryProvider extends AbstractProvider<ModelFactory<PDMPModel>> {

	private Provider<AdaptiveMSHRN> hrnProvider;

	@Inject
	public AdaptiveMSHRNModelFactoryProvider(HierarchicalConfiguration config, Provider<AdaptiveMSHRN> hrnProvider) {
		super(config, "ModelParameters");
		this.hrnProvider = hrnProvider;
	}

	@Override
	public ModelFactory<PDMPModel> get() {
		return new ModelFactory<PDMPModel>() {

			@Override
			public PDMPModel createModel() {
				AdaptiveMSHRN hrn = hrnProvider.get();
				AdaptiveMSHRNModel hrnModel = new AdaptiveMSHRNModel(hrn);
				hrnModel.setExposeOptionalState(config().getBoolean("exposeOptionalState", false));
				return hrnModel;
			}

		};
	}

}
