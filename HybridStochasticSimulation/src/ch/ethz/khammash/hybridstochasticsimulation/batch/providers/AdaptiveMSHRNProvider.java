package ch.ethz.khammash.hybridstochasticsimulation.batch.providers;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.averaging.AveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.networks.AdaptiveMSHRN;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork;

import com.google.inject.Inject;

public class AdaptiveMSHRNProvider extends AbstractProvider<AdaptiveMSHRN> {

	private MSHybridReactionNetwork hrn;
	private AveragingUnit averagingUnit;

	@Inject
	public AdaptiveMSHRNProvider(HierarchicalConfiguration config, MSHybridReactionNetwork hrn, AveragingUnit averagingUnit) {
		super(config, "ModelParameters");
		this.hrn = hrn;
		this.averagingUnit = averagingUnit;
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
		adaptiveHrn.setPrintMessages(config().getBoolean("printMessages", false));
		adaptiveHrn.setAveragingUnit(averagingUnit);
		return adaptiveHrn;
	}

}
