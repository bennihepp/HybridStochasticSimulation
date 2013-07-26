package ch.ethz.khammash.hybridstochasticsimulation.graphs;

import java.util.List;

import org.jgrapht.DirectedGraph;

public interface ReactionNetworkGraph extends DirectedGraph<SpeciesVertex, ReactionEdge> {

    SpeciesVertex getSpeciesVertex(int species);

    List<ReactionEdge> getReactionEdges(int reaction);

}
