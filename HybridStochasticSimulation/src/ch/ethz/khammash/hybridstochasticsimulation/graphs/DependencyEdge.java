package ch.ethz.khammash.hybridstochasticsimulation.graphs;

public class DependencyEdge {

	private SpeciesVertex source;
	private SpeciesVertex target;

	public DependencyEdge(SpeciesVertex source, SpeciesVertex target) {
		this.source = source;
		this.target = target;
	}

	public SpeciesVertex getSource() {
		return source;
	}

	public SpeciesVertex getTarget() {
		return target;
	}

	@Override
	public String toString() {
		return source + "->" + target;
	}

}
