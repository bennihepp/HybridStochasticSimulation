package ch.ethz.khammash.hybridstochasticsimulation.injection.guiceproviders;

import java.io.Serializable;

import javax.inject.Provider;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.models.AdaptiveMSHRNModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.AdaptiveMSHRN;

import com.google.inject.Inject;

public class AdaptiveMSHRNModelProvider extends AbstractProvider<PDMPModel> implements Serializable {

	private static final long serialVersionUID = 1L;

	//	private Provider<AdaptiveMSHRN> hrnProvider;
	private AdaptiveMSHRN baseHrn;

	@Inject
	public AdaptiveMSHRNModelProvider(HierarchicalConfiguration config, Provider<AdaptiveMSHRN> hrnProvider) {
		super(config, "ModelParameters");
//		this.hrnProvider = hrnProvider;
		baseHrn = hrnProvider.get();
	}

	@Override
	public PDMPModel get() {
//		AdaptiveMSHRN hrn = hrnProvider.get();
		AdaptiveMSHRN hrn = AdaptiveMSHRN.createCopy(baseHrn);
		AdaptiveMSHRNModel hrnModel = new AdaptiveMSHRNModel(hrn);
		hrnModel.setExposeOptionalState(config().getBoolean("exposeOptionalState", false));
		return hrnModel;
	}

}
