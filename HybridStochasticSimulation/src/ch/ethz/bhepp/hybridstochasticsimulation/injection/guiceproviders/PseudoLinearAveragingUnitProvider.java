package ch.ethz.bhepp.hybridstochasticsimulation.injection.guiceproviders;

import java.util.Set;

import javax.inject.Provider;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.bhepp.hybridstochasticsimulation.averaging.PseudoLinearAveragingUnit;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.MassActionReactionNetwork;

import com.google.inject.Inject;

public class PseudoLinearAveragingUnitProvider extends AbstractAveragingUnitProvider<PseudoLinearAveragingUnit> {

	@Inject
	public PseudoLinearAveragingUnitProvider(HierarchicalConfiguration config, Provider<MassActionReactionNetwork> networkProvider) {
		super(config, networkProvider);
	}

	@Override
	protected PseudoLinearAveragingUnit getAveragingUnit(MassActionReactionNetwork network, Set<Integer> importantSpecies) {
		PseudoLinearAveragingUnit au = new PseudoLinearAveragingUnit(network, importantSpecies);
		au.stopIfAveragingBecomesInvalid(config().getBoolean("stopIfAveragingBecomesInvalid", true));
		au.warnIfAveragingBecomesInvalid(config().getBoolean("warnIfAveragingBecomesInvalid", true));
		au.performPseudoLinearAveragingOnlyOnce(config().getBoolean("performPseudoLinearAveragingOnlyOnce", true));
		return au;
	}

}
