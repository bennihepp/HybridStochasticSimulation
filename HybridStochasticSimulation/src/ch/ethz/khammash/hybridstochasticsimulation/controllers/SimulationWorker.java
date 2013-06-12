package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import java.util.concurrent.Callable;

import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;

public interface SimulationWorker<T extends ReactionNetworkModel, E extends TrajectoryRecorder<T>> extends Callable<E>, Runnable {

	public E simulate();

	public E call();

	public void run();

}
