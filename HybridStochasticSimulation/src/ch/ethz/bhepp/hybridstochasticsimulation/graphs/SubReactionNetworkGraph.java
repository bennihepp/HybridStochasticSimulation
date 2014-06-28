package ch.ethz.bhepp.hybridstochasticsimulation.graphs;

import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jgrapht.graph.DirectedSubgraph;

import ch.ethz.bhepp.hybridstochasticsimulation.networks.MassActionReactionNetwork;


public class SubReactionNetworkGraph extends DirectedSubgraph<SpeciesVertex, ReactionEdge> implements ReactionNetworkGraph {

	private static final long serialVersionUID = 7151147840597310L;

	private Map<Integer, SpeciesVertex> verticesMap;
	private Map<Integer, List<ReactionEdge>> edgesMap;

	public static SubReactionNetworkGraph createSubReactionNetworkGraph(MassActionReactionNetwork network, Map<Integer,SpeciesVertex> verticesMap, Map<Integer,List<ReactionEdge>> edgesMap) {
		Set<SpeciesVertex> speciesSet = new HashSet<>(verticesMap.values());
		Set<ReactionEdge> reactionSet = new HashSet<>();
		for (List<ReactionEdge> edges : edgesMap.values())
			reactionSet.addAll(edges);
		SubReactionNetworkGraph subNetwork = new SubReactionNetworkGraph(network, speciesSet, reactionSet);
		subNetwork.verticesMap = verticesMap;
		subNetwork.edgesMap = edgesMap;
		return subNetwork;
	}

	private SubReactionNetworkGraph(MassActionReactionNetwork network, Set<SpeciesVertex> speciesSet, Set<ReactionEdge> reactionSet) {
		super(network.getGraph(), speciesSet, reactionSet);
	}

    public SpeciesVertex getSpeciesVertex(int species) {
    	return verticesMap.get(species);
    }

    public List<ReactionEdge> getReactionEdges(int reaction) {
    	return Collections.unmodifiableList(edgesMap.get(reaction));
    }

}
