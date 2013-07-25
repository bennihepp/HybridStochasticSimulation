package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import java.util.Iterator;
import java.util.Set;

import ch.ethz.khammash.hybridstochasticsimulation.graphs.ReactionNetworkGraph;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;

import com.google.common.base.Predicate;
import com.google.common.collect.Iterators;
import com.google.common.collect.Sets;

public class MaximumSizeSubnetworksEnumerator implements SubnetworksEnumerator {

	private ReactionNetworkGraph graph;
	private int maxSize;

	public MaximumSizeSubnetworksEnumerator(ReactionNetworkGraph graph, int maxSize) {
		this.graph = graph;
		this.maxSize = maxSize;
	}

	@Override
	public Iterator<Set<SpeciesVertex>> iterator() {
		Set<SpeciesVertex> allSpecies = graph.vertexSet();
		Set<Set<SpeciesVertex>> speciesPowerset = Sets.powerSet(allSpecies);
		Iterator<Set<SpeciesVertex>> baseIterator = speciesPowerset.iterator();
		return Iterators.filter(baseIterator, new Predicate<Set<SpeciesVertex>>() {

			public boolean apply(Set<SpeciesVertex> input) {
				return input.size() <= maxSize;
			}

		});
//		final Iterator<Set<SpeciesVertex>> baseIterator = speciesPowerset.iterator();
//		return new Iterator<Set<SpeciesVertex>>() {
//
//			private Set<SpeciesVertex> next = null;
//
//			@Override
//			public boolean hasNext() {
//				while (baseIterator.hasNext()) {
//					if (next == null)
//						next = baseIterator.next();
//					if (accept(next))
//						return true;
//				}
//				return false;
//			}
//
//			@Override
//			public Set<SpeciesVertex> next() {
//				if (hasNext())
//					return next;
//				else
//					throw new NoSuchElementException();
//			}
//
//			@Override
//			public void remove() {
//				// TODO Auto-generated method stub
//				
//			}
//			
//		}
	}

}
