package ch.ethz.bhepp.ode.nonstiff;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.nonstiff.RungeKuttaIntegrator;
import org.apache.commons.math3.util.FastMath;

import ch.ethz.bhepp.ode.FixedStepSolver;
import ch.ethz.bhepp.ode.Ode;
import ch.ethz.bhepp.ode.Solver.InitializationException;

public class FixedCommonsMathIntegratorAdapter implements FixedStepSolver {

	private double stepSize;
	private RungeKuttaIntegrator integrator;
	private FirstOrderDifferentialEquations commonsMathOde;
	private Ode baseOde;
	private double t;
	private double[] x;
	private double t1;

	public FixedCommonsMathIntegratorAdapter(double stepSize, RungeKuttaIntegrator integrator) {
		this.stepSize = stepSize;
		this.integrator = integrator;
	}

	@Override
	public void initialize(Ode ode) throws InitializationException {
		this.baseOde = ode;
		this.commonsMathOde = new FirstOrderDifferentialEquations() {
			
			@Override
			public int getDimension() {
				return baseOde.getDimensionOfVectorField();
			}
			
			@Override
			public void computeDerivatives(double t, double[] x, double[] xDot)
					throws MaxCountExceededException, DimensionMismatchException {
				baseOde.computeVectorField(t, x, xDot);
			}
		};
	}

	@Override
	public void prepareStep(double t0, double[] x0, double t1) throws IllegalStateException {
		this.t = t0;
		this.x = x0;
		this.t1 = t1;
	}

	@Override
	public double integrateStep() {
		return integrateStep(stepSize);
	}

	@Override
	public double integrateStep(double step) {
		double tNext = t + FastMath.min(stepSize, t1 - this.t);
		double[] xNext = integrator.singleStep(commonsMathOde, t, x, tNext);
		t = tNext;
		for (int i=0; i < xNext.length; i++)
			x[i] = xNext[i];
		return t;
	}

	@Override
	public double integrateStep(double t, double[] x, double t1, double[] xOut) {
		double step = FastMath.min(stepSize, t1 - t);
		return integrateStep(t, x, xOut, step);
	}

	@Override
	public double integrateStep(double t, double[] x, double[] xOut, double step) {
		double tNext = t + step;
		double[] xNext = integrator.singleStep(commonsMathOde, t, x, tNext);
		for (int i=0; i < xNext.length; i++)
			xOut[i] = xNext[i];
		return tNext;
	}

	@Override
	public void setStepSize(double stepSize) throws UnsupportedOperationException {
		this.stepSize = stepSize;
	}

}
