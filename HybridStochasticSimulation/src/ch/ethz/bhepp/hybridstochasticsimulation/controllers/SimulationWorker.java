package ch.ethz.bhepp.hybridstochasticsimulation.controllers;

import java.util.concurrent.Callable;

import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.TrajectoryRecorder;

public interface SimulationWorker extends Callable<TrajectoryRecorder>, Runnable {

	TrajectoryRecorder simulate();

	TrajectoryRecorder call();

	void run();

}
