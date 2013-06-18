package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;

import org.apache.commons.math3.stat.descriptive.StatisticalSummary;

import ch.ethz.khammash.hybridstochasticsimulation.factories.FiniteTrajectoryRecorderFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.ModelFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.RandomDataGeneratorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.SimulatorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.TrajectoryRecorderFactory;
import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;

public interface SimulationController<T extends ReactionNetworkModel, E extends TrajectoryRecorder<T>> {

	void setExecutorService(ExecutorService executor);

	void setSimulatorFactory(SimulatorFactory<Simulator<T, E>> simulatorFactory);

//	void setTrajectoryRecorderFactory(TrajectoryRecorderFactory<E> trajectoryRecorderFactory);

	void setRandomDataGeneratorFactory(RandomDataGeneratorFactory rdgFactory);

//	E simulateTrajectory(T model, double t0, double[] x0, double t1);

	void simulateTrajectory(T model, E tr, double t0, double[] x0, double t1);

	StatisticalSummary[][] simulateTrajectoryDistribution(
			int runs, ModelFactory<T> modelFactory, final TrajectoryRecorderFactory<E> trFactory,
			double[] tSeries, double[] x0)
			throws InterruptedException, CancellationException, ExecutionException;

	StatisticalSummary[][] simulateTrajectoryDistribution(
			int runs, ModelFactory<T> modelFactory, FiniteTrajectoryRecorderFactory<E> trFactory,
			double[] tSeries, double[] x0)
			throws InterruptedException, CancellationException, ExecutionException;

}
