package ch.ethz.bhepp.hybridstochasticsimulation.providers;

import javax.inject.Inject;

import org.apache.commons.math3.analysis.solvers.BrentSolver;
import org.apache.commons.math3.analysis.solvers.UnivariateSolver;


public class BrentSolverProvider implements ObjProvider<UnivariateSolver> {

	private double relativeAccuracy = 1e-6;
	private double absoluteAccuracy = 1e-6;
	private double functionValueAccuracy = 1e-6;

	@Inject
	public BrentSolverProvider() {
	}

	public BrentSolverProvider(double absoluteAccuracy) {
		setAbsoluteAccuracy(absoluteAccuracy);
	}

	public BrentSolverProvider(double relativeAccuracy, double absoluteAccuracy) {
		setRelativeAccuracy(relativeAccuracy);
		setAbsoluteAccuracy(absoluteAccuracy);
	}

	public BrentSolverProvider(double relativeAccuracy, double absoluteAccuracy, double functionValueAccuracy) {
		setRelativeAccuracy(relativeAccuracy);
		setAbsoluteAccuracy(absoluteAccuracy);
		setFunctionValueAccuracy(functionValueAccuracy);
	}

	public double getRelativeAccuracy() {
		return relativeAccuracy;
	}

	public void setRelativeAccuracy(double relativeAccuracy) {
		this.relativeAccuracy = relativeAccuracy;
	}

	public double getAbsoluteAccuracy() {
		return absoluteAccuracy;
	}

	public void setAbsoluteAccuracy(double absoluteAccuracy) {
		this.absoluteAccuracy = absoluteAccuracy;
	}

	public double getFunctionValueAccuracy() {
		return functionValueAccuracy;
	}

	public void setFunctionValueAccuracy(double functionValueAccuracy) {
		this.functionValueAccuracy = functionValueAccuracy;
	}

	@Override
	public UnivariateSolver get() {
		return new BrentSolver(getRelativeAccuracy(), getAbsoluteAccuracy(), getFunctionValueAccuracy());
	}

}
