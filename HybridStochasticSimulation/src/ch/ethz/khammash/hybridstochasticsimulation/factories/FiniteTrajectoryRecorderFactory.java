package ch.ethz.khammash.hybridstochasticsimulation.factories;

import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;

public interface FiniteTrajectoryRecorderFactory extends TrajectoryRecorderFactory {

	FiniteTrajectoryRecorder createTrajectoryRecorder();

}
