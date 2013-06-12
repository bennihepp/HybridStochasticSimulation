package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;

public interface TrajectoryRecorder<T extends ReactionNetworkModel> extends Trajectory {

	public void setModel(T model);

	public void setInitialState(double t0, double[] x0);

	public void setFinalState(double t1, double[] x1);

	public void handleReactionEvent(int reaction, double t, double[] newX);

}
