package ch.ethz.khammash.hybridstochasticsimulation.factories;

import ch.ethz.khammash.ode.Solver;

public interface SolverFactory {

	Solver createSolver();

}
