package ch.ethz.khammash.hybridstochasticsimulation.sandbox;

import org.jgrapht.EdgeFactory;

public class DependencyEdgeFactory implements EdgeFactory<SpeciesVertex, DependencyEdge> {

	public DependencyEdgeFactory() {
	}

	@Override
	public DependencyEdge createEdge(SpeciesVertex source, SpeciesVertex target) {
		return new DependencyEdge(source, target);
	}

}
