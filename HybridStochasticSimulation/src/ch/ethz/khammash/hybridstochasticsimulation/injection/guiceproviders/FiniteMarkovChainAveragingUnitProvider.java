package ch.ethz.khammash.hybridstochasticsimulation.injection.guiceproviders;

import java.util.Set;

import javax.inject.Provider;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.averaging.FiniteMarkovChainAveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;

import com.google.inject.Inject;

public class FiniteMarkovChainAveragingUnitProvider extends AbstractAveragingUnitProvider<FiniteMarkovChainAveragingUnit> {

	@Inject
	public FiniteMarkovChainAveragingUnitProvider(HierarchicalConfiguration config, Provider<UnaryBinaryReactionNetwork> networkProvider) {
		super(config, networkProvider);
	}

	@Override
	protected FiniteMarkovChainAveragingUnit getAveragingUnit(UnaryBinaryReactionNetwork network, Set<SpeciesVertex> importantSpeciesVertices) {
		return new FiniteMarkovChainAveragingUnit(network, importantSpeciesVertices);
	}

}
