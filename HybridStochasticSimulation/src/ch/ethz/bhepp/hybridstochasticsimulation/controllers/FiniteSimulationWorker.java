package ch.ethz.bhepp.hybridstochasticsimulation.controllers;

import java.util.concurrent.Callable;

import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.TrajectoryRecorder;

public interface FiniteSimulationWorker extends Callable<TrajectoryRecorder>, Runnable {

	FiniteTrajectoryRecorder simulate();

	FiniteTrajectoryRecorder call();

	void run();

}
