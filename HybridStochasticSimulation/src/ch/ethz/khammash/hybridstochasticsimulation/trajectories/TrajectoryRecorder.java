package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;

public interface TrajectoryRecorder<T extends ReactionNetworkModel> extends Trajectory {

	void setSimulator(Simulator<T, ? extends TrajectoryRecorder<T>> simulator);

	void setModel(T model);

	void setInitialState(double t0, double[] x0);

	void setInitialState(double t0, double[] x0, int numOfStates);

	void setFinalState(double t1, double[] x1);

	void reportState(double t, double[] x);

//	void handleReactionEvent(int reaction, double t, double[] newX);

}
