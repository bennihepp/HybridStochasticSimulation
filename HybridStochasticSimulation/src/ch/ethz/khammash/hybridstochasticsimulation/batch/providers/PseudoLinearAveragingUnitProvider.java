package ch.ethz.khammash.hybridstochasticsimulation.batch.providers;

import java.util.Set;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.averaging.PseudoLinearAveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;

import com.google.inject.Inject;

public class PseudoLinearAveragingUnitProvider extends AbstractAveragingUnitProvider<PseudoLinearAveragingUnit> {

	@Inject
	public PseudoLinearAveragingUnitProvider(HierarchicalConfiguration config, UnaryBinaryReactionNetwork network) {
		super(config, network);
	}

	@Override
	protected PseudoLinearAveragingUnit getAveragingUnit(double theta, UnaryBinaryReactionNetwork network, Set<SpeciesVertex> importantSpeciesVertices) {
		PseudoLinearAveragingUnit au = new PseudoLinearAveragingUnit(theta, network, importantSpeciesVertices);
		au.stopIfAveragingBecomesInvalid(config().getBoolean("stopIfAveragingBecomesInvalid", true));
		au.warnIfAveragingBecomesInvalid(config().getBoolean("warnIfAveragingBecomesInvalid", true));
		au.performPseudoLinearAveragingOnlyOnce(config().getBoolean("performPseudoLinearAveragingOnlyOnce", true));
		return au;
	}

}
