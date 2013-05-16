package ch.ethz.khammash.hybridstochasticsimulation;

public interface ReactionHandler {

	public void setInitialState(double t, double[] x);

	public void setFinalState(double t, double[] x);

	public void handleReaction(int reaction, double t, double[] newX);

}
