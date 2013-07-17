package ch.ethz.khammash.hybridstochasticsimulation.batch.providers;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.networks.AdaptiveMSHRN;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork;

import com.google.inject.Inject;

public class AdaptiveMSHRNProvider extends AbstractProvider<AdaptiveMSHRN> {

	private MSHybridReactionNetwork hrn;

	@Inject
	public AdaptiveMSHRNProvider(HierarchicalConfiguration config, MSHybridReactionNetwork hrn) {
		super(config, "ModelParameters");
		this.hrn = hrn;
	}

	@Override
	public AdaptiveMSHRN get() {
		AdaptiveMSHRN adaptiveHrn = AdaptiveMSHRN.createFrom(hrn);
		String deltaKey = "delta";
		if (config().getMaxIndex(deltaKey) >= 0)
			adaptiveHrn.setDelta(config().getDouble(deltaKey));
		String epsilonKey = "epsilon";
		if (config().getMaxIndex(epsilonKey) >= 0)
			adaptiveHrn.setDelta(config().getDouble(epsilonKey));
		String xiKey = "xi";
		if (config().getMaxIndex(xiKey) >= 0)
			adaptiveHrn.setDelta(config().getDouble(xiKey));
		return adaptiveHrn;
	}

}
