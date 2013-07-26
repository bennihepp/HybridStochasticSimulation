package ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders;

import java.util.Set;

import javax.inject.Provider;

import org.apache.commons.configuration.HierarchicalConfiguration;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.averaging.ZeroDeficiencyAveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;
import ch.ethz.khammash.hybridstochasticsimulation.providers.ObjProvider;

import com.google.inject.Inject;

public class ZeroDeficiencyAveragingUnitProvider extends AbstractAveragingUnitProvider<ZeroDeficiencyAveragingUnit> {

	private ObjProvider<RandomDataGenerator> rdgFactory;

	@Inject
	public ZeroDeficiencyAveragingUnitProvider(HierarchicalConfiguration config, Provider<UnaryBinaryReactionNetwork> networkProvider,
			ObjProvider<RandomDataGenerator> rdgFactory) {
		super(config, networkProvider);
		this.rdgFactory = rdgFactory;
	}

	@Override
	protected ZeroDeficiencyAveragingUnit getAveragingUnit(UnaryBinaryReactionNetwork network, Set<SpeciesVertex> importantSpeciesVertices) {
		RandomDataGenerator rdg = rdgFactory.get();
		boolean printMessages = config().getBoolean("printMessages", false);
		return new ZeroDeficiencyAveragingUnit(network, importantSpeciesVertices, rdg, printMessages);
	}

}
