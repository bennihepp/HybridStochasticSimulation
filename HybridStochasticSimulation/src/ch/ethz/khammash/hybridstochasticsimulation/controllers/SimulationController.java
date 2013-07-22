package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import java.util.List;
import java.util.concurrent.ExecutorService;

import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.providers.ObjProvider;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteStatisticalSummaryTrajectory;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;

public interface SimulationController<T extends ReactionNetworkModel> {

	void setExecutorService(ExecutorService executor);

	void setSimulatorProvider(ObjProvider<Simulator<T>> simulatorProvider);

	void setRandomDataGeneratorProvider(ObjProvider<RandomDataGenerator> rdgProvider);

//	E simulateTrajectory(T model, double t0, double[] x0, double t1);

	void simulateTrajectory(T model, TrajectoryRecorder tr, double t0, double[] x0, double t1);

	List<TrajectoryRecorder> simulateTrajectories(int runs,
			ObjProvider<? extends T> modelProvider, ObjProvider<? extends TrajectoryRecorder> trProvider,
			double t0, double[] x0, double t1);

	FiniteStatisticalSummaryTrajectory simulateTrajectoryDistribution(
			int runs, ObjProvider<? extends T> modelProvider, ObjProvider<? extends FiniteTrajectoryRecorder> trProvider,
			double t0, double[] x0, double t1);

}
