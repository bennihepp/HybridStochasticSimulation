package ch.ethz.bhepp.hybridstochasticsimulation.injection.guiceproviders;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import javax.inject.Provider;

import org.apache.commons.configuration.HierarchicalConfiguration;
import org.apache.commons.lang3.ArrayUtils;

import ch.ethz.bhepp.hybridstochasticsimulation.averaging.AllSubnetworkEnumerator;
import ch.ethz.bhepp.hybridstochasticsimulation.averaging.ModularAveragingUnit;
import ch.ethz.bhepp.hybridstochasticsimulation.averaging.SubnetworkEnumerator;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.MassActionReactionNetwork;

import com.google.inject.Inject;

public abstract class AbstractAveragingUnitProvider<T extends ModularAveragingUnit> extends AbstractProvider<T> {

	private Provider<MassActionReactionNetwork> networkProvider;

	@Inject
	public AbstractAveragingUnitProvider(HierarchicalConfiguration config, Provider<MassActionReactionNetwork> networkProvider) {
		super(config, "AveragingParameters");
		this.networkProvider = networkProvider;
	}

	@Override
	public T get() {
		MassActionReactionNetwork network = networkProvider.get();
		String importantSpeciesKey = "importantSpecies";
		int[] importantSpeciesIndices = dataConfig().getIntArray(importantSpeciesKey, new int[0]);
		Set<Integer> importantSpecies = new HashSet<>(Arrays.asList(ArrayUtils.toObject(importantSpeciesIndices)));
		T averagingUnit = getAveragingUnit(network, importantSpecies);
		SubnetworkEnumerator subnetworksEnumerator = new AllSubnetworkEnumerator(network);
		// TODO
//		int maxSize = config().getInt("maxSize", -1);
//		if (maxSize > 0) {
//			MaximumSizeSpeciesSubsetFilter filter = new MaximumSizeSpeciesSubsetFilter(maxSize);
//			subnetworksEnumerator = new FilteredSubnetworkEnumerator(subnetworksEnumerator, filter);
//		}
		averagingUnit.setSubnetworkEnumerator(subnetworksEnumerator);
		return averagingUnit;
	}

	protected abstract T getAveragingUnit(MassActionReactionNetwork network, Set<Integer> importantSpecies);

}
