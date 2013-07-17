package ch.ethz.khammash.hybridstochasticsimulation.sandbox;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.ejml.simple.SimpleMatrix;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.alg.StrongConnectivityInspector;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DirectedSubgraph;

import ch.ethz.khammash.hybridstochasticsimulation.networks.ReactionNetwork;

public class ComplexGraph extends DefaultDirectedGraph<ComplexVertex, ComplexEdge> {

	private static final long serialVersionUID = 7992771587636755157L;
	private ConnectivityInspector<ComplexVertex, ComplexEdge> ci;

	private StrongConnectivityInspector<ComplexVertex, ComplexEdge> si;

	public static ComplexGraph createFromReactionNetwork(ReactionNetwork net) {
		return new ComplexGraph(net);
	}

	protected ComplexGraph(ReactionNetwork net) {
		super(new ComplexEdgeFactory());
		init(net);
	}

	private ComplexEdge createComplexEdge(ComplexVertex source, ComplexVertex target) {
		return new ComplexEdge(source, target);
	}

	private ComplexVertex getComplexVertex(int[] stochiometries) {
		ComplexVertex vertex = new ComplexVertex(stochiometries);
		if (!containsVertex(vertex))
			addVertex(vertex);
		return vertex;
//		if (complexVertexSet.contains(vertex)).containsKey(stochiometries))
//			vertex = complexVertexMap.get(stochiometries);
//		else {
//			vertex = new ComplexVertex(stochiometries);
//			addVertex(vertex);
//			complexVertexMap.put(stochiometries, vertex);
//		}
//		return vertex;
	}

    private void addComplexEdge(ComplexVertex consumptionVertex, ComplexVertex productionVertex) {
    	ComplexEdge edge = createComplexEdge(consumptionVertex, productionVertex);
		addEdge(consumptionVertex, productionVertex, edge);
	}

    final private void init(ReactionNetwork net) {
		for (int r=0; r < net.getNumberOfReactions(); r++) {
			int[] consumptionStochiometries = net.getConsumptionStochiometries(r);
			int[] productionStochiometries = net.getProductionStochiometries(r);
			ComplexVertex consumptionVertex = getComplexVertex(consumptionStochiometries);
			ComplexVertex productionVertex = getComplexVertex(productionStochiometries);
			addComplexEdge(consumptionVertex, productionVertex);
		}
    }

	public int computeDeficiency() {
		int rank;
		if (vertexSet().size() > 0) {
			ComplexVertex cv = vertexSet().iterator().next();
			int numOfSpecies = cv.getStochiometries().length;
			SimpleMatrix matrix = new SimpleMatrix(vertexSet().size(), numOfSpecies);
			int i = 0;
			for (ComplexVertex v : vertexSet()) {
				int[] stochiometries = v.getStochiometries();
				for (int s=0; s < numOfSpecies; s++)
					matrix.set(i, s, stochiometries[s]);
				i++;
			}
			rank = matrix.svd().rank();
		} else
			rank = 0;

		List<Set<ComplexVertex>> connectedSets = getConnectedSets();

		int numOfComplexes = vertexSet().size();
		int dimOfStochiometricSubspace = rank;
		int numOfLinkageClasses = connectedSets.size();
		int deficiency = numOfComplexes - numOfLinkageClasses - dimOfStochiometricSubspace;

		for (Set<ComplexVertex> connectedSet : connectedSets) {
			System.out.println("  Connected set:");
			for (ComplexVertex v : connectedSet) {
				System.out.println("    " + v);
			}
		}

		System.out.println("  Number of complexes: " + numOfComplexes);
		System.out.println("  Dimension of stochiometric subspace: " + dimOfStochiometricSubspace);
		System.out.println("  Number of linkage classes: " + numOfLinkageClasses);
		System.out.println("  Deficiency: " + deficiency);
		System.out.println();

		return deficiency;
	}

	// A chemical reaction network is weakly reversible if and only if all connected components of the complex graph are strongly connected
	public boolean isWeaklyReversible() {
		List<Set<ComplexVertex>> connectedSets = getConnectedSets();
		for (Set<ComplexVertex> connectedSet : connectedSets) {
			Set<ComplexEdge> connectedSetEdges = new HashSet<>();
			for (ComplexEdge e : this.edgeSet()) {
				if (connectedSet.contains(getEdgeSource(e)) && connectedSet.contains(getEdgeTarget(e)))
					connectedSetEdges.add(e);
			}
			DirectedSubgraph<ComplexVertex, ComplexEdge> subgraph = new DirectedSubgraph<>(this, connectedSet, connectedSetEdges);
			StrongConnectivityInspector<ComplexVertex, ComplexEdge> sci = new StrongConnectivityInspector<>(subgraph);
			if (!sci.isStronglyConnected())
				return false;
		}
		return true;
	}

    public List<Set<ComplexVertex>> getConnectedSets() {
    	initCi();
    	return ci.connectedSets();
    }

    public List<Set<ComplexVertex>> getStronglyConnectedSets() {
    	initSi();
    	return si.stronglyConnectedSets();
    }

    private void initCi() {
    	if (ci == null)
    		ci = new ConnectivityInspector<>(this);
    }

    private void initSi() {
    	if (si == null)
    		si = new StrongConnectivityInspector<>(this);
    }

}
