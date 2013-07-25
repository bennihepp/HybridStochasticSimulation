package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import java.util.Iterator;
import java.util.Set;

import ch.ethz.khammash.hybridstochasticsimulation.graphs.ReactionNetworkGraph;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;

import com.google.common.collect.Sets;

public class AllSubnetworksEnumerator implements SubnetworksEnumerator {

	private ReactionNetworkGraph graph;

	public AllSubnetworksEnumerator(ReactionNetworkGraph graph) {
		this.graph = graph;
	}

	@Override
	public Iterator<Set<SpeciesVertex>> iterator() {
		Set<SpeciesVertex> allSpecies = graph.vertexSet();
		Set<Set<SpeciesVertex>> speciesPowerset = Sets.powerSet(allSpecies);
		return speciesPowerset.iterator();
	}

}
