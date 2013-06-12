package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;

public class DefaultSimulationWorker<T extends ReactionNetworkModel, E extends TrajectoryRecorder<T>>
		implements SimulationWorker<T, E> {

	private Simulator<T, E> simulator;
	private T model;
	private E tr;
	private double t0;
	private double[] x0;
	private double t1;

	public DefaultSimulationWorker(
			Simulator<T, E> simulator,
			T model,
			E tr,
			double t0, double[] x0, double t1) {
		this.simulator = simulator;
		this.model = model;
		this.tr = tr;
		this.t0 = t0;
		this.x0 = x0;
		this.t1 = t1;
	}

	@Override
	public E simulate() {
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
	public E call() {
		return simulate();
	}

}
