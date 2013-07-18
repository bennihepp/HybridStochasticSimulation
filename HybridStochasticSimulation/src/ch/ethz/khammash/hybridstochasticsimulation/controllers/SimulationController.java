package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import java.util.List;
import java.util.concurrent.ExecutorService;

import ch.ethz.khammash.hybridstochasticsimulation.factories.FiniteTrajectoryRecorderFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.ModelFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.RandomDataGeneratorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.SimulatorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.TrajectoryRecorderFactory;
import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteStatisticalSummaryTrajectory;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;

public interface SimulationController<T extends ReactionNetworkModel> {

	void setExecutorService(ExecutorService executor);

	void setSimulatorFactory(SimulatorFactory<T> simulatorFactory);

//	void setTrajectoryRecorderFactory(TrajectoryRecorderFactory<E> trajectoryRecorderFactory);

	void setRandomDataGeneratorFactory(RandomDataGeneratorFactory rdgFactory);

//	E simulateTrajectory(T model, double t0, double[] x0, double t1);

	void simulateTrajectory(T model, TrajectoryRecorder tr, double t0, double[] x0, double t1);

	List<TrajectoryRecorder> simulateTrajectories(int runs,
			ModelFactory<T> modelFactory, TrajectoryRecorderFactory trFactory,
			double t0, double[] x0, double t1);

	FiniteStatisticalSummaryTrajectory simulateTrajectoryDistribution(
			int runs, ModelFactory<T> modelFactory, FiniteTrajectoryRecorderFactory trFactory,
			double t0, double[] x0, double t1);

}
