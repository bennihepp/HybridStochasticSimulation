package ch.ethz.khammash.hybridstochasticsimulation.batch.providers;

import java.util.HashSet;
import java.util.Set;

import javax.inject.Provider;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.averaging.AveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;

import com.google.inject.Inject;

public abstract class AbstractAveragingUnitProvider<T extends AveragingUnit> extends AbstractProvider<T> {

	private Provider<UnaryBinaryReactionNetwork> networkProvider;

	@Inject
	public AbstractAveragingUnitProvider(HierarchicalConfiguration config, Provider<UnaryBinaryReactionNetwork> networkProvider) {
		super(config, "AveragingParameters");
		this.networkProvider = networkProvider;
	}

	@Override
	public T get() {
		UnaryBinaryReactionNetwork network = networkProvider.get();
		double theta = config().getDouble("theta");
		int[] importantSpecies = dataConfig().getIntArray("importantSpecies");
		Set<SpeciesVertex> importantSpeciesVertices = new HashSet<>();
		for (int s : importantSpecies)
			importantSpeciesVertices.add(network.getGraph().getSpeciesVertex(s));
		return getAveragingUnit(theta, network, importantSpeciesVertices);
	}

	protected abstract T getAveragingUnit(double theta, UnaryBinaryReactionNetwork network,
			Set<SpeciesVertex> importantSpeciesVertices);

}
