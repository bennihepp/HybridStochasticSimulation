package ch.ethz.khammash.hybridstochasticsimulation.factories;

import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;

public interface TrajectoryRecorderFactory<E extends TrajectoryRecorder<?>> {

	E createTrajectoryRecorder();

}
