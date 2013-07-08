package ch.ethz.khammash.hybridstochasticsimulation.factories;

import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPSimulator;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;
import ch.ethz.khammash.ode.Solver;

public class PDMPSimulatorFactory
		implements SimulatorFactory<Simulator<PDMPModel>> {

	private SolverFactory solverFactory;

	public PDMPSimulatorFactory(SolverFactory solverFactory) {
		this.solverFactory = solverFactory;
	}

	public void setSolverFactory(SolverFactory solverFactory) {
		this.solverFactory = solverFactory;
	}

	@Override
	public PDMPSimulator createSimulator(RandomDataGenerator rdg) {
		Solver solver = solverFactory.createSolver();
		PDMPSimulator sim = new PDMPSimulator(solver, rdg);
		return sim;
	}

}