package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import java.util.Iterator;
import java.util.Set;

import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;

import com.google.common.collect.Sets;

public class AllSubnetworkEnumerator implements SubnetworkEnumerator {

	private UnaryBinaryReactionNetwork network;
//	private Predicate<SubnetworkDescription> subnetworkFilter;

	public AllSubnetworkEnumerator(UnaryBinaryReactionNetwork network) {
		this.network = network;
//		this(network, Predicates.<SubnetworkDescription>alwaysTrue());
	}

//	public AllSubnetworksEnumerator(UnaryBinaryReactionNetwork network, Predicate<SubnetworkDescription> subnetworkFilter) {
//		this.network = network;
//		this.subnetworkFilter = subnetworkFilter;
//	}

	@Override
	public Iterator<SubnetworkDescription> iterator() {
		Set<SpeciesVertex> allSpecies = network.getGraph().vertexSet();
		Set<Set<SpeciesVertex>> speciesPowerset = Sets.powerSet(allSpecies);
		final Iterator<Set<SpeciesVertex>> baseIterator = speciesPowerset.iterator();
		Iterator<SubnetworkDescription> iterator = new Iterator<SubnetworkDescription>() {

			@Override
			public boolean hasNext() {
				return baseIterator.hasNext();
			}

			@Override
			public SubnetworkDescription next() {
				return new SubnetworkDescription(baseIterator.next(), network);
			}

			@Override
			public void remove() {
				// TODO Auto-generated method stub
				
			}
		};
//		return Iterators.filter(iterator, subnetworkFilter);
		return iterator;
	}

}
