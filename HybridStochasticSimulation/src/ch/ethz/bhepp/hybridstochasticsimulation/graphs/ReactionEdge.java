package ch.ethz.bhepp.hybridstochasticsimulation.graphs;

public class ReactionEdge {

	private SpeciesVertex source;
	private SpeciesVertex target;
	private int reaction;
	private double kappa;

	public ReactionEdge(SpeciesVertex source, SpeciesVertex target) {
		this(-1, Double.NaN, source, target);
	}

	public ReactionEdge(int reaction, SpeciesVertex source, SpeciesVertex target) {
		this(reaction, Double.NaN, source, target);
	}

	public ReactionEdge(int reaction, double kappa, SpeciesVertex source, SpeciesVertex target) {
		this.reaction = reaction;
		this.setKappa(kappa);
		this.source = source;
		this.target = target;
	}

	public int getReaction() {
		return reaction;
	}

	public void setReaction(int reaction) {
		this.reaction = reaction;
	}

	public double getKappa() {
		return kappa;
	}

	public void setKappa(double kappa) {
		this.kappa = kappa;
	}

	public SpeciesVertex getSource() {
		return source;
	}

	public SpeciesVertex getTarget() {
		return target;
	}

	@Override
	public String toString() {
		return source + "->" + target + " [" + reaction + "]";
	}

}
