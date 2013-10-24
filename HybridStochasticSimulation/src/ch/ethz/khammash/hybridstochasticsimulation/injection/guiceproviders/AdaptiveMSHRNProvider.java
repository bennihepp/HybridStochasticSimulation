package ch.ethz.khammash.hybridstochasticsimulation.injection.guiceproviders;

import javax.inject.Provider;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.networks.AdaptiveMSHRN;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork;

import com.google.inject.Inject;

public class AdaptiveMSHRNProvider extends AbstractProvider<AdaptiveMSHRN> {

	private Provider<MSHybridReactionNetwork> hrnProvider;

	@Inject
	public AdaptiveMSHRNProvider(HierarchicalConfiguration config, Provider<MSHybridReactionNetwork> hrnProvider) {
		super(config, "ModelParameters");
		this.hrnProvider = hrnProvider;
	}

	@Override
	public AdaptiveMSHRN get() {
		MSHybridReactionNetwork hrn = hrnProvider.get();
		AdaptiveMSHRN adaptiveHrn = AdaptiveMSHRN.createFrom(hrn);
		String deltaKey = "delta";
		if (config().getMaxIndex(deltaKey) >= 0)
			adaptiveHrn.setDelta(config().getDouble(deltaKey));
		String epsilonKey = "eta";
		if (config().getMaxIndex(epsilonKey) >= 0)
			adaptiveHrn.setEta(config().getDouble(epsilonKey));
		String xiKey = "xi";
		if (config().getMaxIndex(xiKey) >= 0)
			adaptiveHrn.setXi(config().getDouble(xiKey));
		String thetaKey = "theta";
		if (config().getMaxIndex(thetaKey) >= 0)
			adaptiveHrn.setTheta(config().getDouble(thetaKey));
		adaptiveHrn.setLogMessages(config().getBoolean("printMessages", false));
		return adaptiveHrn;
	}

}
