package ch.ethz.bhepp.ode.nonstiff;

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import ch.ethz.bhepp.ode.EventFunction;
import ch.ethz.bhepp.ode.EventObserver;
import ch.ethz.bhepp.ode.FiniteTimepointProvider;
import ch.ethz.bhepp.ode.Ode;
import ch.ethz.bhepp.ode.StateObserver;

public class RungeKutta4thOrderSolver extends EulerSolver {

	protected DenseMatrix64F xV;
	protected DenseMatrix64F kV;
	protected DenseMatrix64F tmpV;
	protected DenseMatrix64F dummyV;
	protected double[][] k;

	public RungeKutta4thOrderSolver(double step) {
		super(step);
	}

	@Override
	public void initialize(Ode ode, EventFunction ef, StateObserver stateObserver, EventObserver eventObserver) {
		super.initialize(ode, ef, stateObserver, eventObserver);
		xV = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		kV = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		tmpV = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		dummyV = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		this.k = new double[c.length][ode.getDimensionOfVectorField()];
	}

	@Override
	protected double findNextState(double t, double[] x, double step, double maxStep) {
		if (step > maxStep)
			step = maxStep;
		findNextStateExplicitFixed(t, x, step, x);
		return t + step;
	}

	private final double[][] a = { { }, { 1/2. }, { 0., 1/2. }, { 0., 0., 1. } };
	private final double[] b = { 1/6., 2/6., 2/6., 1/6. };
	private final double[] c = { 0., 1/2., 1/2., 1. };

	protected void findNextStateExplicitFixed(double t, double[] x, double step, double[] out) {
//		ode.computeVectorField(t, x, k[0]);
//		xV.setData(x);
//		for (int i=1; i < k.length; i++) {
//			tmpV.set(xV);
//			for (int j=0; j < i; j++) {
//				kV.setData(k[j]);
//				CommonOps.add(tmpV, step * a[i][j], kV, tmpV);
//			}
//			ode.computeVectorField(t + step * c[i], tmpV.getData(), k[i]);
//		}
//		dummyV.setData(out);
//		dummyV.set(xV);
//		for (int i=0; i < k.length; i++) {
//			kV.setData(k[i]);
//			CommonOps.add(dummyV, step * b[i], kV, dummyV);
//		}
		// Manual
		xV.setData(x);
		tmpV.set(xV);
		//
		ode.computeVectorField(t, x, k[0]);
		kV.setData(k[0]);
		CommonOps.scale(step, kV);
		//
		kV.setData(k[0]);
		CommonOps.add(tmpV, a[1][0], kV, tmpV);
		ode.computeVectorField(t + step * c[1], tmpV.getData(), k[1]);
		kV.setData(k[1]);
		CommonOps.scale(step, kV);
		//
		tmpV.set(xV);
		kV.setData(k[0]);
		CommonOps.add(tmpV, a[2][0], kV, tmpV);
		kV.setData(k[1]);
		CommonOps.add(tmpV, a[2][1], kV, tmpV);
		ode.computeVectorField(t + step * c[2], tmpV.getData(), k[2]);
		kV.setData(k[2]);
		CommonOps.scale(step, kV);
		//
		tmpV.set(xV);
		kV.setData(k[0]);
		CommonOps.add(tmpV, a[3][0], kV, tmpV);
		kV.setData(k[1]);
		CommonOps.add(tmpV, a[3][1], kV, tmpV);
		kV.setData(k[2]);
		CommonOps.add(tmpV, a[3][2], kV, tmpV);
		ode.computeVectorField(t + step * c[3], tmpV.getData(), k[3]);
		kV.setData(k[3]);
		CommonOps.scale(step, kV);
		//
		dummyV.setData(out);
		dummyV.set(xV);
		kV.setData(k[0]);
		CommonOps.add(dummyV, b[0], kV, dummyV);
		kV.setData(k[1]);
		CommonOps.add(dummyV, b[1], kV, dummyV);
		kV.setData(k[2]);
		CommonOps.add(dummyV, b[2], kV, dummyV);
		kV.setData(k[3]);
		CommonOps.add(dummyV, b[3], kV, dummyV);
	}

//	@Override
//	protected void findNextStateExplicitFixed(double t, double[] x, double step, double[] out) {
//		ode.computeVectorField(t, x, k1);
//		for (int i=0; i < k1.length; i++)
//			tmp[i] = x[i] + step / 2.0 * k1[i];
//		ode.computeVectorField(t, tmp, k2);
//		for (int i=0; i < k1.length; i++)
//			tmp[i] = x[i] + step / 2.0 * k2[i];
//		ode.computeVectorField(t, tmp, k3);
//		for (int i=0; i < k1.length; i++)
//			tmp[i] = x[i] + step * k3[i];
//		ode.computeVectorField(t, tmp, k4);
//		for (int i=0; i < k1.length; i++)
//			out[i] = x[i] + 1 / 6.0 * step * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
//	}

    public static void main(String args[]) {
    	MyOde ode = new MyOde();
    	MyEf ef = new MyEf();
    	MyObserver observer = new MyObserver();
    	MyEventObserver eventObserver = new MyEventObserver();
    	double step = 0.1;
    	RungeKutta4thOrderSolver q = new RungeKutta4thOrderSolver(step);
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
