package ch.ethz.bhepp.hybridstochasticsimulation.injection.guiceproviders;

import java.util.Set;

import javax.inject.Provider;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.bhepp.hybridstochasticsimulation.averaging.FiniteMarkovChainAveragingUnit;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.MassActionReactionNetwork;

import com.google.inject.Inject;

public class FiniteMarkovChainAveragingUnitProvider extends AbstractAveragingUnitProvider<FiniteMarkovChainAveragingUnit> {

	@Inject
	public FiniteMarkovChainAveragingUnitProvider(HierarchicalConfiguration config, Provider<MassActionReactionNetwork> networkProvider) {
		super(config, networkProvider);
	}

	@Override
	protected FiniteMarkovChainAveragingUnit getAveragingUnit(MassActionReactionNetwork network, Set<Integer> importantSpecies) {
		return new FiniteMarkovChainAveragingUnit(network, importantSpecies);
	}

}
