package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import java.util.Iterator;
import java.util.Set;

import ch.ethz.khammash.hybridstochasticsimulation.graphs.ReactionNetworkGraph;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;

import com.google.common.base.Predicate;
import com.google.common.base.Predicates;
import com.google.common.collect.Iterators;
import com.google.common.collect.Sets;

public class FilteredSubnetworksEnumerator implements SubnetworkEnumerator {

	private ReactionNetworkGraph graph;
	private Predicate<Set<SpeciesVertex>> subnetworkFilter;

	public FilteredSubnetworksEnumerator(ReactionNetworkGraph graph) {
		this(graph, Predicates.<Set<SpeciesVertex>>alwaysTrue());
	}

	public FilteredSubnetworksEnumerator(ReactionNetworkGraph graph, Predicate<Set<SpeciesVertex>> subnetworkFilter) {
		this.graph = graph;
		this.subnetworkFilter = subnetworkFilter;
	}

	@Override
	public Iterator<Set<SpeciesVertex>> iterator() {
		Set<SpeciesVertex> allSpecies = graph.vertexSet();
		Set<Set<SpeciesVertex>> speciesPowerset = Sets.powerSet(allSpecies);
		Iterator<Set<SpeciesVertex>> baseIterator = speciesPowerset.iterator();
		return Iterators.filter(baseIterator, subnetworkFilter);
	}

}
