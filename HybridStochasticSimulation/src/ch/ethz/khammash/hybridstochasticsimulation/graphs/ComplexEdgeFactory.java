package ch.ethz.khammash.hybridstochasticsimulation.graphs;

import org.jgrapht.EdgeFactory;

public class ComplexEdgeFactory implements EdgeFactory<ComplexVertex, ComplexEdge> {

	public ComplexEdgeFactory() {
	}

	@Override
	public ComplexEdge createEdge(ComplexVertex source, ComplexVertex target) {
		return new ComplexEdge(source, target);
	}

}
