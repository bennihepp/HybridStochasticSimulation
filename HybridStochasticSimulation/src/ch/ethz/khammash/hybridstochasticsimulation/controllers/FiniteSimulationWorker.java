package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import java.util.concurrent.Callable;

import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;

public interface FiniteSimulationWorker extends Callable<TrajectoryRecorder>, Runnable {

	FiniteTrajectoryRecorder simulate();

	FiniteTrajectoryRecorder call();

	void run();

}
