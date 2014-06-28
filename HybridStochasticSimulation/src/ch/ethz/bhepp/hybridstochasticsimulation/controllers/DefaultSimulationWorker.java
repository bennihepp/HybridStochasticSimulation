package ch.ethz.bhepp.hybridstochasticsimulation.controllers;

import ch.ethz.bhepp.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.bhepp.hybridstochasticsimulation.simulators.Simulator;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.TrajectoryRecorder;

public class DefaultSimulationWorker<T extends ReactionNetworkModel>
	implements SimulationWorker {

	private Simulator<T> simulator;
	private T model;
	private TrajectoryRecorder tr;
	private double t0;
	private double[] x0;
	private double t1;

	public DefaultSimulationWorker(
			Simulator<T> simulator,
			T model,
			TrajectoryRecorder tr,
			double t0, double[] x0, double t1) {
		this.simulator = simulator;
		this.model = model;
		this.tr = tr;
		this.t0 = t0;
		this.x0 = x0;
		this.t1 = t1;
	}

	@Override
	public TrajectoryRecorder simulate() {
		double[] x1 = new double[x0.length];
		simulator.addTrajectoryRecorder(tr);
		simulator.simulate(model, t0, x0, t1, x1);
		simulator.clearTrajectoryRecorders();
		return tr;
	}

	@Override
	public void run() {
		simulate();
	}

	@Override
	public TrajectoryRecorder call() {
		return simulate();
	}

}
