package ch.ethz.khammash.hybridstochasticsimulation.sandbox;

import java.util.List;

import org.jgrapht.graph.DefaultDirectedGraph;


public class DependencyGraph extends DefaultDirectedGraph<SpeciesVertex, DependencyEdge> {

	private static final long serialVersionUID = 4205342458526263526L;

	private ReactionNetworkGraph graph;

	public DependencyGraph(ReactionNetworkGraph graph) {
		super(new DependencyEdgeFactory());
		this.graph = graph;
		init();
	}

    final private void init() {
    	for (SpeciesVertex v : graph.vertexSet()) {
    		addVertex(v);
    	}
    	for (SpeciesVertex source : vertexSet()) {
    		List<SpeciesVertex> vertices = graph.getReachableVertices(source);
    		for (SpeciesVertex target : vertices)
    			addEdge(source, target);
    	}
    }

    public ReactionNetworkGraph getReactionNetworkGraph() {
    	return graph;
    }

}
