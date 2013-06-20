package ch.ethz.khammash.hybridstochasticsimulation.simulators;

import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;

public interface Simulator<T extends ReactionNetworkModel, E extends TrajectoryRecorder<T>> {

	double simulate(T model, final double t0, final double[] x0, double t1, double[] x1);

	void addTrajectoryRecorder(E tr);

	void removeTrajectoryRecorder(E tr);

	void clearTrajectoryRecorders();

}
