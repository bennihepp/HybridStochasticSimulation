package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;

public interface FiniteTrajectoryRecorderFactory<E extends TrajectoryRecorder<?>> {

	public E createTrajectoryRecorder(double[] tSeries);

}
