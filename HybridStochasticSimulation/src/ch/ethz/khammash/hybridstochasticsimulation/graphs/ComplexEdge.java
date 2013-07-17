package ch.ethz.khammash.hybridstochasticsimulation.graphs;

public class ComplexEdge {

	private ComplexVertex source;
	private ComplexVertex target;

	public ComplexEdge(ComplexVertex source, ComplexVertex target) {
		this.source = source;
		this.target = target;
	}

	public ComplexVertex getSource() {
		return source;
	}

	public ComplexVertex getTarget() {
		return target;
	}

	@Override
	public String toString() {
		return source + "->" + target;
	}

}
