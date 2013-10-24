package ch.ethz.khammash.hybridstochasticsimulation.graphs;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.jgrapht.alg.StrongConnectivityInspector;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DirectedSubgraph;
import org.jgrapht.traverse.BreadthFirstIterator;

import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;


public class DefaultReactionNetworkGraph extends DefaultDirectedGraph<SpeciesVertex, ReactionEdge> implements ReactionNetworkGraph {

	private static final long serialVersionUID = -8643146341228206681L;

	private ArrayList<SpeciesVertex> vertices;
	private ArrayList<LinkedList<ReactionEdge>> edgeLists;
	private StrongConnectivityInspector<SpeciesVertex, ReactionEdge> si;

	public DefaultReactionNetworkGraph(UnaryBinaryReactionNetwork network) {
		super(new ReactionEdgeFactory());
		init(network);
	}

	private ReactionEdge createReactionEdge(int reaction, String label, double kappa, SpeciesVertex source, SpeciesVertex target) {
		return new ReactionEdge(reaction, kappa, source, target);
	}

    final private void init(UnaryBinaryReactionNetwork network) {
		vertices = new ArrayList<SpeciesVertex>(network.getNumberOfSpecies());
		edgeLists = new ArrayList<LinkedList<ReactionEdge>>(network.getNumberOfReactions());
		for (int s=0; s < network.getNumberOfSpecies(); s++) {
			SpeciesVertex v = new SpeciesVertex(s, network.getSpeciesLabel(s));
			vertices.add(v);
    		addVertex(v);
		}

    	for (int r=0; r < network.getNumberOfReactions(); r++) {
    		LinkedList<ReactionEdge> edgeList = new LinkedList<ReactionEdge>();
    		edgeLists.add(edgeList);
    		for (int s1=0; s1 < network.getNumberOfSpecies(); s1++) {
    			if (network.getConsumptionStoichiometry(s1, r) > 0) {
        			SpeciesVertex v1 = vertices.get(s1);
    	    		for (int s2=0; s2 < network.getNumberOfSpecies(); s2++) {
            			SpeciesVertex v2 = vertices.get(s2);
    	    			if (network.getProductionStoichiometry(s2, r) > 0) {
    	    				ReactionEdge edge = createReactionEdge(r, network.getReactionLabel(r), network.getRateParameter(r), v1, v2);
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
