package ch.ethz.bhepp.hybridstochasticsimulation.simulators;

import ch.ethz.bhepp.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.TrajectoryRecorder;

public interface Simulator<T extends ReactionNetworkModel> {

	double simulate(T model, final double t0, final double[] x0, double t1, double[] x1);

	void addTrajectoryRecorder(TrajectoryRecorder tr);

	void removeTrajectoryRecorder(TrajectoryRecorder tr);

	void clearTrajectoryRecorders();

	void setPrintMessages(boolean printMessages);

	void setShowProgress(boolean showProgress);

}
