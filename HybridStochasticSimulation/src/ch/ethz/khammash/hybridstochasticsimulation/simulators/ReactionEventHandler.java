package ch.ethz.khammash.hybridstochasticsimulation.simulators;

public interface ReactionEventHandler {

	public void setInitialState(double t, double[] x);

	public void setFinalState(double t, double[] x);

	public void handleReactionEvent(int reaction, double t, double[] newX);

}
