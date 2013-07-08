package ch.ethz.khammash.hybridstochasticsimulation.factories;

import ch.ethz.khammash.ode.Solver;
import ch.ethz.khammash.ode.cvode.CVodeSolver;

public class CVodeSolverFactory implements SolverFactory {

	private double relTolerance;
	private double absTolerance;

	public CVodeSolverFactory() {
		this(1e-3, 1e-3);
	}

	public CVodeSolverFactory(double relTolerance, double absTolerance) {
		this.setRelTolerance(relTolerance);
		this.setAbsTolerance(absTolerance);
	}

	@Override
	public Solver createSolver() {
		return new CVodeSolver(getRelTolerance(), getAbsTolerance());
	}

	public double getRelTolerance() {
		return relTolerance;
	}

	public void setRelTolerance(double relTolerance) {
		this.relTolerance = relTolerance;
	}

	public double getAbsTolerance() {
		return absTolerance;
	}

	public void setAbsTolerance(double absTolerance) {
		this.absTolerance = absTolerance;
	}
}