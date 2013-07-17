package ch.ethz.khammash.hybridstochasticsimulation.sandbox;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.jgrapht.alg.StrongConnectivityInspector;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DirectedSubgraph;
import org.jgrapht.traverse.BreadthFirstIterator;

import ch.ethz.khammash.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;


public class ReactionNetworkGraph extends DefaultDirectedGraph<SpeciesVertex, ReactionEdge> {

	private static final long serialVersionUID = -8643146341228206681L;

	private ArrayList<SpeciesVertex> vertices;
	private ArrayList<LinkedList<ReactionEdge>> edgeLists;
	private StrongConnectivityInspector<SpeciesVertex, ReactionEdge> si;

	public ReactionNetworkGraph(DefaultUnaryBinaryReactionNetwork net) {
		super(new ReactionEdgeFactory());
		init(net);
	}

    final private void init(DefaultUnaryBinaryReactionNetwork net) {
		vertices = new ArrayList<SpeciesVertex>(net.getNumberOfSpecies());
		edgeLists = new ArrayList<LinkedList<ReactionEdge>>(net.getNumberOfReactions());
		for (int s=0; s < net.getNumberOfSpecies(); s++) {
			SpeciesVertex v = new SpeciesVertex(s);
			vertices.add(v);
    		addVertex(v);
		}

    	for (int r=0; r < net.getNumberOfReactions(); r++) {
    		LinkedList<ReactionEdge> edgeList = new LinkedList<ReactionEdge>();
    		edgeLists.add(edgeList);
    		for (int s1=0; s1 < net.getNumberOfSpecies(); s1++) {
    			if (net.getConsumptionStochiometry(s1, r) > 0) {
        			SpeciesVertex v1 = vertices.get(s1);
    	    		for (int s2=0; s2 < net.getNumberOfSpecies(); s2++) {
            			SpeciesVertex v2 = vertices.get(s2);
    	    			if (net.getProductionStochiometry(s2, r) > 0) {
    	    				ReactionEdge edge = new ReactionEdge(r, net.getRateParameter(r), v1, v2);
    	    				edgeList.add(edge);
    	    				addEdge(v1, v2, edge);
    	    			}
    	    		}
    			}
    		}
    	}
    }

    public List<SpeciesVertex> getReachableVertices(SpeciesVertex source) {
    	BreadthFirstIterator<SpeciesVertex, ReactionEdge> it = new BreadthFirstIterator<SpeciesVertex, ReactionEdge>(this, source);
    	LinkedList<SpeciesVertex> reachableVertices = new LinkedList<SpeciesVertex>();
    	while (it.hasNext()) {
    		SpeciesVertex v = it.next();
    		reachableVertices.add(v);
    	}
    	return new ArrayList<SpeciesVertex>(reachableVertices);
    }

    public SpeciesVertex getSpeciesVertex(int species) {
    	return vertices.get(species);
    }

    public List<ReactionEdge> getReactionEdges(int reaction) {
    	return Collections.unmodifiableList(edgeLists.get(reaction));
    }

    public List<Set<SpeciesVertex>> getStronglyConnectedSets() {
    	initSi();
    	return si.stronglyConnectedSets();
    }

    public List<DirectedSubgraph<SpeciesVertex, ReactionEdge>> getStronglyConnectedSubgraphs() {
    	initSi();
    	return si.stronglyConnectedSubgraphs();
    }

    private void initSi() {
    	if (si == null)
    		si = new StrongConnectivityInspector<SpeciesVertex, ReactionEdge>(this);
    }

}
