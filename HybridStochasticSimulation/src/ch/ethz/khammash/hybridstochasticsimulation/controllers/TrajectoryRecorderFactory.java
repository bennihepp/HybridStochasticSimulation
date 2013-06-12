package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;

public interface TrajectoryRecorderFactory<E extends TrajectoryRecorder<?>> {

	public E createTrajectoryRecorder();

}
