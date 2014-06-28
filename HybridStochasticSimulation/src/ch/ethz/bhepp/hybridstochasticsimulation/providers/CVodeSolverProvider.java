package ch.ethz.bhepp.hybridstochasticsimulation.providers;

import javax.inject.Inject;

import ch.ethz.bhepp.ode.Solver;
import ch.ethz.bhepp.ode.cvode.CVodeSolver;

public class CVodeSolverProvider implements ObjProvider<Solver> {

	private double relTolerance;
	private double absTolerance;

	@Inject
	public CVodeSolverProvider() {
		this(1e-3, 1e-3);
	}

	public CVodeSolverProvider(double relTolerance, double absTolerance) {
		this.setRelTolerance(relTolerance);
		this.setAbsTolerance(absTolerance);
	}

	@Override
	public Solver get() {
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