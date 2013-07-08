package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import java.util.concurrent.Callable;

import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;

public interface SimulationWorker extends Callable<TrajectoryRecorder>, Runnable {

	TrajectoryRecorder simulate();

	TrajectoryRecorder call();

	void run();

}
