package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;

public class DefaultFiniteSimulationWorker<T extends ReactionNetworkModel> implements FiniteSimulationWorker {

	private DefaultSimulationWorker<T> worker;

	public DefaultFiniteSimulationWorker(
			Simulator<T> simulator,
			T model,
			FiniteTrajectoryRecorder tr,
			double t0, double[] x0, double t1) {
		worker = new DefaultSimulationWorker<T>(simulator, model, tr, t0, x0, t1);
	}

	@Override
	public FiniteTrajectoryRecorder simulate() {
		return (FiniteTrajectoryRecorder)worker.simulate();
	}

	@Override
	public void run() {
		worker.run();
	}

	@Override
	public FiniteTrajectoryRecorder call() {
		return simulate();
	}

}
