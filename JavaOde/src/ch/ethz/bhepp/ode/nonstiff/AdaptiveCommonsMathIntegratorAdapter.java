package ch.ethz.bhepp.ode.nonstiff;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.nonstiff.RungeKuttaIntegrator;
import org.apache.commons.math3.util.FastMath;
import org.ejml.ops.CommonOps;

import ch.ethz.bhepp.ode.AdaptiveStepSolver;
import ch.ethz.bhepp.ode.Ode;
import ch.ethz.bhepp.ode.Solver.InitializationException;

public class AdaptiveCommonsMathIntegratorAdapter implements AdaptiveStepSolver {

	private double absTolerance;
	private double relTolerance;
	private double stepSize;
	private RungeKuttaIntegrator integrator;
	private FirstOrderDifferentialEquations commonsMathOde;
	private Ode baseOde;
	private double t;
	private double[] x;
	private double t1;
	private double[] xPrev;
	private double tPrev;

	public AdaptiveCommonsMathIntegratorAdapter(double stepSize, RungeKuttaIntegrator integrator) {
		absTolerance = 1e-3;
		relTolerance = 1e-3;
		this.stepSize = stepSize;
		this.integrator = integrator;
	}

	public double getAbsTolerance() {
		return absTolerance;
	}

	public void setAbsTolerance(double absTolerance) {
		this.absTolerance = absTolerance;
	}

	public double getRelTolerance() {
		return relTolerance;
	}

	public void setRelTolerance(double relTolerance) {
		this.relTolerance = relTolerance;
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
		this.tPrev = t0;
		this.xPrev = new double[commonsMathOde.getDimension()];
		this.t1 = t1;
	}

	@Override
	public void prepareStep(double t0, double[] x0, double t1, double dt0) throws IllegalStateException {
		setCurrentStepSize(dt0);
		prepareStep(t0, x0, t1);
	}

	@Override
	public double integrateStep() {
		return integrateStep(x);
	}

	private final double TINY = 1e-30;

	@Override
	public double integrateStep(double[] xOut) {
		tPrev = t;
		for (int i=0; i < xPrev.length; i++)
			xPrev[i] = x[i];
		double curStep = stepSize;
		double newStep = stepSize;
		double[] xTemp1;
		double[] xTemp2;
		while (true) {
			curStep = newStep;
			curStep = FastMath.min(curStep, t1 - this.t);
			double halfStep = curStep / 2.0;
			xTemp1 = integrator.singleStep(commonsMathOde, t, x, t + curStep);
			xTemp2 = integrator.singleStep(commonsMathOde, t, x, t + halfStep);
			xTemp2 = integrator.singleStep(commonsMathOde, t + halfStep, xTemp2, t + curStep);
			double err = Double.NaN;
			if (absTolerance > 0) {
				double absError = 0.0;
				for (int i=0; i < commonsMathOde.getDimension(); i++)
					absError = FastMath.max(FastMath.abs(xTemp2[i] - xTemp1[i]), absError);
				err = absError / absTolerance;
			}
			if (relTolerance > 0) {
				double relError = 0.0;
				for (int i=0; i < commonsMathOde.getDimension(); i++) {
					double val = FastMath.abs((xTemp2[i] - xTemp1[i]) / x[i]) + TINY;
					if (!Double.isNaN(val))
						relError = FastMath.max(val, relError);
				}
				relError = relError / relTolerance;
				if (relError > err || Double.isNaN(err))
					err = relError;
			}
			if (Double.isNaN(err) || Double.isInfinite(err)) {
				newStep = curStep / 2.0;
				continue;
			}
			if (err == 0.0) {
				newStep = 2 * curStep;
				break;
			}
			if (err <= 1.0) {
				newStep = FastMath.min(0.9 * curStep / err, 5.0 * curStep);
				break;
			} else
				newStep = FastMath.max(0.9 * curStep / err, 0.1 * curStep);
				newStep = 0.9 * curStep / err;
		}
		stepSize = newStep;
		t = t + curStep;
		for (int i=0; i < xTemp2.length; i++)
			x[i] = xTemp2[i];
		return t;
	}

	@Override
	public void computeInterpolatedSolution(double t, double[] xOut) throws IllegalStateException {
		double deltaT = this.t - tPrev;
		double factor = (t - tPrev) / deltaT;
		for (int i=0; i < xPrev.length; i++)
			xOut[i] = x[i] + factor * (x[i] - xPrev[i]);
	}

	@Override
	public double getCurrentStepSize() throws IllegalStateException, UnsupportedOperationException {
		return stepSize;
	}

	@Override
	public double getLastStepSize() throws IllegalStateException, UnsupportedOperationException {
		throw new UnsupportedOperationException();
	}

	@Override
	public void setCurrentStepSize(double stepSize) throws IllegalStateException, UnsupportedOperationException {
		this.stepSize = stepSize;
	}

	@Override
	public void setMinimumStepSize(double minStepSize) throws IllegalStateException, UnsupportedOperationException {
		throw new UnsupportedOperationException();
	}

	@Override
	public void setMaximumStepSize(double maxStepSize) throws IllegalStateException, UnsupportedOperationException {
		throw new UnsupportedOperationException();
	}

}
