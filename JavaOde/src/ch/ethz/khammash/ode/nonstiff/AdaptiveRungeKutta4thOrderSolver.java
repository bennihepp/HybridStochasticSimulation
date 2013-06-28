package ch.ethz.khammash.ode.nonstiff;

import org.apache.commons.math3.util.FastMath;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import ch.ethz.khammash.ode.EventFunction;
import ch.ethz.khammash.ode.EventObserver;
import ch.ethz.khammash.ode.Ode;
import ch.ethz.khammash.ode.StateObserver;

public class AdaptiveRungeKutta4thOrderSolver extends EulerSolver {

	private double relTolerance;
	private double absTolerance;
	private DenseMatrix64F xV;
	private DenseMatrix64F kV;
	private DenseMatrix64F out1V;
	private DenseMatrix64F out2V;
	private DenseMatrix64F errV;
	private DenseMatrix64F tmpV;
	protected DenseMatrix64F dummyV;
	protected double[][] k;


	public AdaptiveRungeKutta4thOrderSolver(double relTolerance, double absTolerance) {
		this(0.1, relTolerance, absTolerance);
	}

	public AdaptiveRungeKutta4thOrderSolver(double step, double relTolerance, double absTolerance) {
		super(step);
		this.relTolerance = relTolerance;
		this.absTolerance = absTolerance;
	}

	@Override
	public void initialize(Ode ode, EventFunction ef, StateObserver stateObserver, EventObserver eventObserver) {
		super.initialize(ode, ef, stateObserver, eventObserver);
		xV = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		kV = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		out1V = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		out2V = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		errV = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		tmpV = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		dummyV = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		this.k = new double[c.length][ode.getDimensionOfVectorField()];
	}

	private final double[][] a = { { }, { 1/4. }, { 3/32., 9/32. }, { 1932/2197., -7200/2197., 7296/2197. },
			{ 439/216., -8., 3680/513., -845/4104. }, { -8/27., 2., -3544/2565., 1859/4104., -11/40. } };
	private final double[] b = { 16/135., 0., 6656/12825., 28561/56430., -9/50., 2/55. };
	private final double[] b_star = { 25/216., 0., 1408/2565., 2197/4104., -1/5., 0. };
	private final double[] c = { 0., 1/4., 3/8., 12/13., 1., 1/2. };

	protected void findNextStateExplicitFixed(double t, double[] x, double step, double[] out1, double[] out2) {
		ode.computeVectorField(t, x, k[0]);
		xV.setData(x);
		for (int i=1; i < k.length; i++) {
			tmpV.set(xV);
			for (int j=0; j < i; j++) {
				kV.setData(k[j]);
				CommonOps.add(tmpV, step * a[i][j], kV, tmpV);
			}
			ode.computeVectorField(t + step * c[i], tmpV.getData(), k[i]);
		}
		dummyV.setData(out1);
		dummyV.set(xV);
		for (int i=0; i < k.length; i++) {
			kV.setData(k[i]);
			CommonOps.add(dummyV, step * b_star[i], kV, dummyV);
		}
		dummyV.setData(out2);
		dummyV.set(xV);
		for (int i=0; i < k.length; i++) {
			kV.setData(k[i]);
			CommonOps.add(dummyV, step * b[i], kV, dummyV);
		}
	}

	private final double TINY = 1e-30;

	protected double findNextState(double t, double[] x, double step, double maxStep) {
		xV.setData(x);
		double curStep = step;
		double newStep = step;
		while (true) {
			curStep = newStep;
//			System.out.println("step: " + step);
//			double halfStep = curStep / 2.0;
			findNextStateExplicitFixed(t, x, curStep, out1V.getData(), out2V.getData());
//			findNextStateExplicitFixed(t, x, halfStep, out2V.getData());
//			findNextStateExplicitFixed(t, out2V.getData(), halfStep, out2V.getData());
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
				if (Double.isInfinite(relErr))
					relErr = TINY;
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
		};
//		if (Double.isNaN(step))
//			throw new Exception("Can't integrate ODE: Step size is NaN");
		if (curStep > maxStep) {
			curStep = maxStep;
			findNextStateExplicitFixed(t, x, curStep, out2V.getData());
		} else
			CommonOps.add(out2V, errV, xV);
		this.step = newStep;
//		System.out.println("done: " + step);
		return t + curStep;
	}

    public static void main(String args[]) {
    	MyOde ode = new MyOde();
    	MyEf ef = new MyEf();
    	MyObserver observer = new MyObserver();
    	MyEventObserver eventObserver = new MyEventObserver();
    	double step = 0.1;
    	AdaptiveRungeKutta4thOrderSolver q = new AdaptiveRungeKutta4thOrderSolver(step, 1e-10, 1e-10);
    	q.initialize(ode, ef, observer, eventObserver);
    	double[] x0 = { 0 };
    	double t0 = 0.0;
    	double t1 = 1.0;
    	q.integrate(t0, x0, t1);
    	x0[0] = 0.0;
    	double[] tSeries = { 0, 0.5, 1.0 };
    	Timepoints timepointProvider = new Timepoints(tSeries);
    	q.integrate(timepointProvider, x0);
        q.dispose();
    };

}
