package ch.ethz.bhepp.hybridstochasticsimulation.injection.guiceproviders;

import javax.inject.Provider;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.bhepp.hybridstochasticsimulation.models.PDMPMSHRNModel;
import ch.ethz.bhepp.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.MSHybridReactionNetwork;

import com.google.inject.Inject;

public class MSHybridReactionNetworkModelProvider extends AbstractObjProvider<PDMPModel> {

	private Provider<MSHybridReactionNetwork> hrnProvider;

	@Inject
	public MSHybridReactionNetworkModelProvider(HierarchicalConfiguration config, Provider<MSHybridReactionNetwork> hrnProvider) {
		super(config, "ModelParameters");
		this.hrnProvider = hrnProvider;
	}

	@Override
	public PDMPModel get() {
		MSHybridReactionNetwork hrn = hrnProvider.get();
		PDMPMSHRNModel hrnModel = new PDMPMSHRNModel(hrn);
		hrnModel.setExposeOptionalState(config().getBoolean("exposeOptionalState", false));
		return hrnModel;
	}

}
