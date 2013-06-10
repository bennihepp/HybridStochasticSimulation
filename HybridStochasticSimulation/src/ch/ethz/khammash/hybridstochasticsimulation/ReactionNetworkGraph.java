package ch.ethz.khammash.hybridstochasticsimulation;

import java.util.List;
import java.util.Set;


import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.alg.StrongConnectivityInspector;

import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetwork;



public class ReactionNetworkGraph extends DefaultDirectedGraph<Integer, DefaultEdge> {

	private static final long serialVersionUID = -3348628197042070775L;

	protected Integer[] vertices;

	public ReactionNetworkGraph(ReactionNetwork net) {
		super(DefaultEdge.class);
		init(net);
	}

    void init(ReactionNetwork net) {
		vertices = new Integer[net.getNumberOfSpecies()];
		for (int s=0; s < vertices.length; s++) {
			vertices[s] = s;
    		addVertex(vertices[s]);
		}

    	for (int r=0; r < net.getNumberOfReactions(); r++)
    		for (int s1=0; s1 < net.getNumberOfSpecies(); s1++)
    			if (net.getConsumptionStochiometry(s1, r) > 0)
    	    		for (int s2=0; s2 < net.getNumberOfSpecies(); s2++)
    	    			if (net.getProductionStochiometry(s2, r) > 0)
    	    				addEdge(vertices[s1], vertices[s2]);
    }

    public List<Set<Integer>> stronglyConnectedSets() {
		StrongConnectivityInspector<Integer, DefaultEdge> si = new StrongConnectivityInspector<Integer, DefaultEdge>(this);
    	return si.stronglyConnectedSets();
    }

}
