package ch.ethz.bhepp.ode.nonstiff;

import org.apache.commons.math3.util.FastMath;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import ch.ethz.bhepp.ode.EventFunction;
import ch.ethz.bhepp.ode.EventObserver;
import ch.ethz.bhepp.ode.FiniteTimepointProvider;
import ch.ethz.bhepp.ode.Ode;
import ch.ethz.bhepp.ode.StateObserver;

public class AdaptiveEulerSolver extends EulerSolver {

	private double relTolerance;
	private double absTolerance;
	private double minStep;
	private DenseMatrix64F xV;
	private DenseMatrix64F out1V;
	private DenseMatrix64F out2V;
	private DenseMatrix64F errV;
	private DenseMatrix64F tmpV;

	public AdaptiveEulerSolver(double relTolerance, double absTolerance) {
		this(0.1, relTolerance, absTolerance, 0.0);
	}

	public AdaptiveEulerSolver(double step, double relTolerance, double absTolerance) {
		this(step, relTolerance, absTolerance, 0.0);
	}

	public AdaptiveEulerSolver(double step, double relTolerance, double absTolerance, double minStep) {
		super(step);
		this.relTolerance = relTolerance;
		this.absTolerance = absTolerance;
		this.minStep = minStep;
	}

	@Override
	public void initialize(Ode ode, EventFunction ef, StateObserver stateObserver, EventObserver eventObserver) {
		super.initialize(ode, ef, stateObserver, eventObserver);
		xV = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		out1V = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		out2V = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		errV = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		tmpV = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
	}

	private final double TINY = 1e-30;

	protected double findNextState(double t, double[] x, double step, double maxStep) {
		xV.setData(x);
		double curStep = step;
		double newStep = step;
		while (true) {
			curStep = newStep;
			double halfStep = curStep / 2.0;
			findNextStateExplicitFixed(t, x, curStep, out1V.getData());
			findNextStateExplicitFixed(t, x, halfStep, out2V.getData());
			findNextStateExplicitFixed(t + halfStep, out2V.getData(), halfStep, out2V.getData());
			CommonOps.sub(out2V, out1V, errV);
//			double errNorm = NormOps.fastNormP2(errV);
			double err = Double.NaN;
			if (absTolerance > 0) {
				double absErr = CommonOps.elementMaxAbs(errV);
				err = absErr / absTolerance;
			}
			if (relTolerance > 0) {
				CommonOps.elementDiv(errV, xV, tmpV);
				double relErr = CommonOps.elementMaxAbs(tmpV) + TINY;
				relErr = relErr / relTolerance;
				if (relErr > err || Double.isNaN(err))
					err = relErr;
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
			if (newStep < minStep)
				break;
		};
//		if (Double.isNaN(step))
//			throw new Exception("Can't integrate ODE: Step size is NaN");
		if (curStep > maxStep) {
			curStep = maxStep;
			findNextStateExplicitFixed(t, x, curStep, out2V.getData());
		} else if (curStep < minStep) {
			curStep = minStep;
			findNextStateExplicitFixed(t, x, curStep, out2V.getData());
		}
//		else
//			CommonOps.add(out2V, errV, xV);
		xV.set(out2V);
		this.stepSize = newStep;
		return t + curStep;
	}

    public static void main(String args[]) {
    	MyOde ode = new MyOde();
    	MyEf ef = new MyEf();
    	MyObserver observer = new MyObserver();
    	MyEventObserver eventObserver = new MyEventObserver();
    	double step = 0.1;
    	AdaptiveEulerSolver q = new AdaptiveEulerSolver(step, 1e-2, 1e-2);
    	q.initialize(ode, ef, observer, eventObserver);
    	double[] x0 = { 0 };
    	double t0 = 0.0;
    	double t1 = 1.0;
    	q.integrate(t0, x0, t1);
    	x0[0] = 0.0;
    	double[] tSeries = { 0, 0.5, 1.0 };
    	FiniteTimepointProvider timepointProvider = new FiniteTimepointProvider(tSeries);
    	q.integrate(timepointProvider, x0);
        q.dispose();
    };

}
