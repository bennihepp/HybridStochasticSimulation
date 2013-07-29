package ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.inject.Provider;

import org.apache.commons.configuration.HierarchicalConfiguration;
import org.apache.commons.lang3.ArrayUtils;

import ch.ethz.khammash.hybridstochasticsimulation.averaging.AveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.averaging.MaximumSizeSubnetworkEnumerator;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;
import ch.ethz.khammash.hybridstochasticsimulation.networks.ReactionNetworkUtils;
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
		int[] importantSpeciesIndices = new int[0];
		String importantSpeciesKey = "importantSpecies";
		String[] importantSpecies = null;
		if (config().getMaxIndex(importantSpeciesKey) >= 0)
				importantSpecies = config().getStringArray(importantSpeciesKey);
		if (importantSpecies != null) {
			List<Integer> importantSpeciesList = new ArrayList<>(importantSpecies.length);
			for (String name : importantSpecies) {
				importantSpeciesList.add(ReactionNetworkUtils.getSpeciesIndex(network, name));
			}
			importantSpeciesIndices = ArrayUtils.toPrimitive(importantSpeciesList.toArray(new Integer[0]));
		} else {
			importantSpeciesIndices = dataConfig().getIntArray("importantSpeciesIndices", importantSpeciesIndices);
		}
		Set<SpeciesVertex> importantSpeciesVertices = new HashSet<>();
		for (int s : importantSpeciesIndices)
			importantSpeciesVertices.add(network.getGraph().getSpeciesVertex(s));
		T averagingUnit = getAveragingUnit(network, importantSpeciesVertices);
		int maxSize = config().getInt("maxSize", -1);
		if (maxSize > 0) {
			MaximumSizeSubnetworkEnumerator subnetworksEnumerator = new MaximumSizeSubnetworkEnumerator(network.getGraph(), maxSize);
			averagingUnit.setSubnetworkEnumerator(subnetworksEnumerator);
		}
		return averagingUnit;
	}

	protected abstract T getAveragingUnit(UnaryBinaryReactionNetwork network,
			Set<SpeciesVertex> importantSpeciesVertices);

}
