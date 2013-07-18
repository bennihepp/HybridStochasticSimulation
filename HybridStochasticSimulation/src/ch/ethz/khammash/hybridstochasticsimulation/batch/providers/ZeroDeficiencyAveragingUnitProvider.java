package ch.ethz.khammash.hybridstochasticsimulation.batch.providers;

import java.util.Set;

import org.apache.commons.configuration.HierarchicalConfiguration;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.averaging.ZeroDeficiencyAveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.factories.RandomDataGeneratorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;

import com.google.inject.Inject;

public class ZeroDeficiencyAveragingUnitProvider extends AbstractAveragingUnitProvider<ZeroDeficiencyAveragingUnit> {

	private RandomDataGeneratorFactory rdgFactory;

	@Inject
	public ZeroDeficiencyAveragingUnitProvider(HierarchicalConfiguration config, UnaryBinaryReactionNetwork network,
			RandomDataGeneratorFactory rdgFactory) {
		super(config, network);
		this.rdgFactory = rdgFactory;
	}

	@Override
	protected ZeroDeficiencyAveragingUnit getAveragingUnit(double theta,
			UnaryBinaryReactionNetwork network,
			Set<SpeciesVertex> importantSpeciesVertices) {
		RandomDataGenerator rdg = rdgFactory.createRandomDataGenerator();
		boolean printMessages = config().getBoolean("printMessages", false);
		return new ZeroDeficiencyAveragingUnit(theta, network, importantSpeciesVertices, rdg, printMessages);
	}

}
