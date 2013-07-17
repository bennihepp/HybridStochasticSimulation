package ch.ethz.khammash.hybridstochasticsimulation.sandbox;

import org.jgrapht.EdgeFactory;

public class ReactionEdgeFactory implements EdgeFactory<SpeciesVertex, ReactionEdge> {

	public ReactionEdgeFactory() {
	}

	@Override
	public ReactionEdge createEdge(SpeciesVertex source, SpeciesVertex target) {
		return new ReactionEdge(source, target);
	}

}
