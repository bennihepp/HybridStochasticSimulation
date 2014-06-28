package ch.ethz.bhepp.hybridstochasticsimulation.graphs;

import java.util.List;
import java.util.Set;

import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.alg.StrongConnectivityInspector;
import org.jgrapht.graph.DefaultDirectedGraph;

import ch.ethz.bhepp.hybridstochasticsimulation.networks.MassActionReactionNetwork;

public class ComplexGraph extends DefaultDirectedGraph<ComplexVertex, ComplexEdge> {

	private static final long serialVersionUID = 7992771587636755157L;
	private ConnectivityInspector<ComplexVertex, ComplexEdge> ci;

	private StrongConnectivityInspector<ComplexVertex, ComplexEdge> si;

	public static ComplexGraph createFromReactionNetwork(MassActionReactionNetwork net) {
		return new ComplexGraph(net);
	}

	protected ComplexGraph(MassActionReactionNetwork net) {
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
	}

    private void addComplexEdge(ComplexVertex consumptionVertex, ComplexVertex productionVertex) {
    	ComplexEdge edge = createComplexEdge(consumptionVertex, productionVertex);
		addEdge(consumptionVertex, productionVertex, edge);
	}

    final private void init(MassActionReactionNetwork net) {
		for (int r=0; r < net.getNumberOfReactions(); r++) {
			int[] consumptionStochiometries = net.getConsumptionStoichiometries(r);
			int[] productionStochiometries = net.getProductionStoichiometries(r);
			ComplexVertex consumptionVertex = getComplexVertex(consumptionStochiometries);
			ComplexVertex productionVertex = getComplexVertex(productionStochiometries);
			addComplexEdge(consumptionVertex, productionVertex);
		}
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
