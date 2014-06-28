package ch.ethz.bhepp.hybridstochasticsimulation.controllers;

import java.util.List;

import ch.ethz.bhepp.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.bhepp.hybridstochasticsimulation.providers.ObjProvider;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FiniteStatisticalSummaryTrajectory;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FiniteTrajectoryMapper;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.TrajectoryRecorder;

public interface SimulationController<T extends ReactionNetworkModel> {

//	void setExecutorService(ExecutorService executor);

//	void setSimulatorProvider(ObjProvider<? extends Simulator<T>> simulatorProvider);

	void simulateTrajectory(T model, TrajectoryRecorder tr, double t0, double[] x0, double t1);

	List<TrajectoryRecorder> simulateTrajectories(int runs,
			ObjProvider<? extends T> modelProvider, ObjProvider<? extends TrajectoryRecorder> trProvider,
			double t0, double[] x0, double t1)
					throws InterruptedException;

	FiniteStatisticalSummaryTrajectory simulateTrajectoryDistribution(
			int runs, ObjProvider<? extends T> modelProvider, ObjProvider<? extends FiniteTrajectoryRecorder> trProvider,
			double t0, double[] x0, double t1)
					throws InterruptedException;

	FiniteStatisticalSummaryTrajectory simulateTrajectoryDistribution(
			int runs, ObjProvider<? extends T> modelProvider, ObjProvider<? extends FiniteTrajectoryRecorder> trProvider,
			FiniteTrajectoryMapper mapper, double t0, double[] x0, double t1)
					throws InterruptedException;

}
