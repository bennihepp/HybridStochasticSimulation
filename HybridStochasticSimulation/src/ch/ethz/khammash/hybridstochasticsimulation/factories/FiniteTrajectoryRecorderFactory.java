package ch.ethz.khammash.hybridstochasticsimulation.factories;

import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;

public interface FiniteTrajectoryRecorderFactory<E extends TrajectoryRecorder<?>> {

	E createTrajectoryRecorder(double[] tSeries);

}
