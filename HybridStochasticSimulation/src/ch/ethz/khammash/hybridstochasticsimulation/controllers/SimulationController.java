package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;

import org.apache.commons.math3.stat.descriptive.StatisticalSummary;

import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;

public interface SimulationController<T extends ReactionNetworkModel, E extends TrajectoryRecorder<T>> {

	public void setExecutorService(ExecutorService executor);

	public void setSimulatorFactory(SimulatorFactory<Simulator<T, E>> simulatorFactory);

//	public void setTrajectoryRecorderFactory(TrajectoryRecorderFactory<E> trajectoryRecorderFactory);

	public void setRandomDataGeneratorFactory(RandomDataGeneratorFactory rdgFactory);

//	public E simulateTrajectory(T model, double t0, double[] x0, double t1);

	public void simulateTrajectory(T model, E tr, double t0, double[] x0, double t1);

	public StatisticalSummary[][] simulateTrajectoryDistribution(
			int runs, ModelFactory<T> modelFactory, FiniteTrajectoryRecorderFactory<E> trFactory,
			double[] tSeries, double[] x0)
			throws InterruptedException, CancellationException, ExecutionException;

}
