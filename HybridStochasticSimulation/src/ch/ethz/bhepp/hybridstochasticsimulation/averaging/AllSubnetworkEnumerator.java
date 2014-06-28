package ch.ethz.bhepp.hybridstochasticsimulation.averaging;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import ch.ethz.bhepp.hybridstochasticsimulation.networks.MassActionReactionNetwork;

import com.google.common.collect.Sets;

public class AllSubnetworkEnumerator implements SubnetworkEnumerator {

	private MassActionReactionNetwork network;
//	private Predicate<SubnetworkDescription> subnetworkFilter;

	public AllSubnetworkEnumerator(MassActionReactionNetwork network) {
		this.network = network;
//		this(network, Predicates.<SubnetworkDescription>alwaysTrue());
	}

//	public AllSubnetworksEnumerator(UnaryBinaryReactionNetwork network, Predicate<SubnetworkDescription> subnetworkFilter) {
//		this.network = network;
//		this.subnetworkFilter = subnetworkFilter;
//	}

	@Override
	public Iterator<SubnetworkDescription> iterator() {
		Set<Integer> allReactions = new HashSet<>(network.getNumberOfReactions());
		for (int r=0; r < network.getNumberOfReactions(); r++)
			allReactions.add(r);
		Set<Set<Integer>> reactionsPowerset = Sets.powerSet(allReactions);
		final Iterator<Set<Integer>> baseIterator = reactionsPowerset.iterator();
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
