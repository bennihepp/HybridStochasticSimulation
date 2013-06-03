package ch.ethz.khammash.hybridstochasticsimulation.simulators;

public interface ReactionEventHandler {

	public void setInitialState(double t0, double[] x0);

	public void setFinalState(double t1, double[] x1);

	public void handleReactionEvent(int reaction, double t, double[] newX);

}
