package ch.ethz.khammash.hybridstochasticsimulation.sandbox;

public class ReactionEdge {

	private SpeciesVertex source;
	private SpeciesVertex target;
	private int reaction;
	private double kappa;

	public ReactionEdge(SpeciesVertex v1, SpeciesVertex v2) {
		this(-1, Double.NaN, v1, v2);
	}

	public ReactionEdge(int reaction, SpeciesVertex v1, SpeciesVertex v2) {
		this(reaction, Double.NaN, v1, v2);
	}

	public ReactionEdge(int reaction, double kappa, SpeciesVertex source, SpeciesVertex target) {
		this.reaction = reaction;
		this.source = source;
		this.target = target;
		this.setKappa(kappa);
	}

	public SpeciesVertex getSource() {
		return source;
	}

	public SpeciesVertex getTarget() {
		return target;
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

	@Override
	public String toString() {
		return source + "->" + target + " [" + reaction + "]";
	}

}
