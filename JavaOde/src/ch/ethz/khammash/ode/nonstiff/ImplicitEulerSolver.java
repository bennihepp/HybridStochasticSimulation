package ch.ethz.khammash.ode.nonstiff;

import org.apache.commons.math3.util.FastMath;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import org.ejml.ops.NormOps;
import org.ejml.ops.SpecializedOps;

import ch.ethz.khammash.ode.EventFunction;
import ch.ethz.khammash.ode.EventObserver;
import ch.ethz.khammash.ode.Ode;
import ch.ethz.khammash.ode.StateObserver;

public class ImplicitEulerSolver extends EulerSolver {

	private double relTolerance;
	private double absTolerance;
	private double rootTolerance;
//	private double[] xDotTmp;
	private DenseMatrix64F xV;
	private DenseMatrix64F xDotV;
	private DenseMatrix64F out1V;
	private DenseMatrix64F out2V;
	private DenseMatrix64F errV;
	private DenseMatrix64F y1;
	private DenseMatrix64F y2;
	private DenseMatrix64F dy;
	private DenseMatrix64F f1;
	private DenseMatrix64F f2;
	private DenseMatrix64F df;
	private DenseMatrix64F tmpV;
	private DenseMatrix64F dummy;
	private DenseMatrix64F B;
	private DenseMatrix64F B_inv;
	private DenseMatrix64F tmpM;

	public ImplicitEulerSolver(double step, double relTolerance, double absTolerance, double rootTolerance) {
		super(step);
		this.relTolerance = relTolerance;
		this.absTolerance = absTolerance;
		this.rootTolerance = rootTolerance;
	}

	@Override
	public void initialize(Ode ode, EventFunction ef,
			StateObserver stateObserver, EventObserver eventObserver) {
		super.initialize(ode, ef, stateObserver, eventObserver);
//		xDotTmp = new double[ode.getDimensionOfVectorField()];
		xV = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		xDotV = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		out1V = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		out2V = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		errV = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		y1 = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		y2 = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		dy = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		f1 = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		f2 = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		df = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		tmpV = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		dummy = new DenseMatrix64F(ode.getDimensionOfVectorField(), 1);
		B = new DenseMatrix64F(ode.getDimensionOfVectorField(), ode.getDimensionOfVectorField());
		B_inv = new DenseMatrix64F(ode.getDimensionOfVectorField(), ode.getDimensionOfVectorField());
		tmpM = new DenseMatrix64F(ode.getDimensionOfVectorField(), ode.getDimensionOfVectorField());
	}

//	protected void findNextState(double t, double[] x, double step) {
//		Array2DRowRealMatrix B = new Array2DRowRealMatrix(x.length, x.length);
//		for (int i=0; i < B.getRowDimension(); i++)
//			for (int j=0; j < B.getColumnDimension(); j++)
//				B.setEntry(i, j, (i == j) ? 1.0e-5 : 0.0);
//		ArrayRealVector v1 = new ArrayRealVector(x);
//		ArrayRealVector f1 = new ArrayRealVector(ode.getDimensionOfVectorField());
//		ArrayRealVector f2 = new ArrayRealVector(ode.getDimensionOfVectorField());
//		ode.computeVectorField(t, x, f1.getDataRef());
//		while (f1.getNorm() > rootTolerance) {
//			RealVector p = B.operate(f1);
//			for (int i=0; i < v1.getDimension(); i++)
//				v1.setEntry(i, v1.getEntry(i) - p.getEntry(i));
//			ode.computeVectorField(t, v1.getDataRef(), f2.getDataRef());
//			ArrayRealVector y = f2.subtract(f1);
//			RealVector oyp = B.operate(y).subtract(p);
//			RealVector pB = B.transpose().operate(p);
//			RealMatrix M = oyp.outerProduct(pB);
//			RealMatrix Q = p.outerProduct(B.operate(y));
//			for (int i=0; i < B.getRowDimension(); i++)
//				for (int j=0; j < B.getColumnDimension(); j++)
//					B.setEntry(i, j, B.getEntry(i, j) - M.getEntry(i, j) / Q.getEntry(i, j));
//			ArrayRealVector fTmp = f1;
//			f1 = f2;
//			f2 = fTmp;
//		}
//		for (int i=0; i < v1.getDimension(); i++)
//			x[i] = v1.getEntry(i);
//	}

	protected void rootFunction(double t, DenseMatrix64F x, DenseMatrix64F y, double step, DenseMatrix64F result) {
//		rootFunctionEuler(t, x, y, step, result);
		rootFunctionTrapezoidal(t, x, y, step, result);
	}

	protected void rootFunctionEuler(double t, DenseMatrix64F x, DenseMatrix64F y, double step, DenseMatrix64F result) {
		ode.computeVectorField(t, y.getData(), xDotV.getData());
		CommonOps.add(x, step, xDotV, result);
		CommonOps.sub(result, y, result);
	}

	protected void rootFunctionTrapezoidal(double t, DenseMatrix64F x, DenseMatrix64F y, double step, DenseMatrix64F result) {
		ode.computeVectorField(t, x.getData(), xDotV.getData());
		ode.computeVectorField(t, y.getData(), result.getData());
		CommonOps.add(xDotV, result, result);
		CommonOps.add(x, 1 / 2.0 * step, result, result);
	}

	protected double findNextState(double t, double[] x, double step, double maxStep) {
//		return findNextStateImplicitVariable(t, x, step, maxStep);
//		return findNextStateImplicitVariableTest(t, x, step, maxStep);
		if (step > maxStep)
			step = maxStep;
		findNextStateImplicitFixed(t, x, step, x);
		return t + step;
	}

	protected double findNextStateImplicitVariable(double t, double[] x, double step, double maxStep) {
		double newStep = step;
		while (true) {
//			System.out.println("step: " + step);
			double halfStep = step / 2.0;
			findNextStateImplicitFixed(t, x, step, out1V.getData());
			findNextStateImplicitFixed(t, x, halfStep, out2V.getData());
			findNextStateImplicitFixed(t, out2V.getData(), halfStep, out2V.getData());
			CommonOps.sub(out2V, out1V, errV);
//			double errNorm = NormOps.fastNormP2(errV);
			double err = Double.NaN;
			if (absTolerance > 0) {
				double absErr = CommonOps.elementMaxAbs(errV);
				err = absErr / absTolerance;
			}
			if (relTolerance > 0) {
				xV.setData(x);
				CommonOps.elementDiv(errV, xV, tmpV);
				double relErr = CommonOps.elementMaxAbs(tmpV);
				relErr = relErr / relTolerance;
				if (relErr > err || Double.isNaN(err))
					err = relErr;
			}
			if (Double.isNaN(err) || Double.isInfinite(err)) {
				step = step / 2.0;
				continue;
			}
			if (err == 0.0) {
				step = 2 * step;
				break;
			}
			if (err <= 1.0) {
				newStep = FastMath.min(0.9 * step / err, 5.0 * step);
				break;
			} else
				newStep = FastMath.max(0.9 * step / err, 0.1 * step);
				newStep = 0.9 * step / err;
		};
//		if (Double.isNaN(step))
//			throw new Exception("Can't integrate ODE: Step size is NaN");
		if (step > maxStep)
			step = maxStep;
		dummy.setData(x);
//		findNextStateImplicitFixed(t, x, prevStep, dummy.getData());
//		CommonOps.add(dummy, errV, dummy);
		CommonOps.add(out2V, errV, dummy);
//		System.arraycopy(out2V.getData(), 0, x, 0, ode.getDimensionOfVectorField());
		this.step = newStep;
//		System.out.println("done: " + step);
		return t + step;
	}

	protected double findNextStateImplicitVariableTest(double t, double[] x, double step, double maxStep) {
		double prevStep;
		while (true) {
//			System.out.println("step: " + step);
			prevStep = step;
			double halfStep = step / 2.0;
			findNextStateImplicitFixed(t, x, step, out1V.getData());
			findNextStateImplicitFixed(t, x, halfStep, out2V.getData());
			findNextStateImplicitFixed(t, out2V.getData(), halfStep, out2V.getData());
			CommonOps.sub(out1V, out2V, errV);
			double errNorm = NormOps.fastNormP2(errV);
			double errRate = errNorm / step;
			if (Double.isNaN(errRate) || Double.isInfinite(errRate)) {
				step = step / 2.0;
				continue;
			}
			if (errRate == 0.0) {
				step = 2 * step;
				break;
			}
			step = 0.9 * step * relTolerance / errRate;
			if (errRate <= relTolerance)
				break;
		};
		if (prevStep > maxStep)
			prevStep = maxStep;
		dummy.setData(x);
		CommonOps.sub(out2V, errV, dummy);
		this.step = step;
		System.out.println("done: " + step);
		return t + prevStep;
	}

	protected void findNextStateImplicitFixed(double t, double[] x, double step, double[] out) {
		CommonOps.setIdentity(B);
		xV.setData(x);
//		System.arraycopy(x, 0, y1.getData(), 0, ode.getDimensionOfVectorField());
		y1.set(xV);
		rootFunctionEuler(t, xV, y1, step, f1);
		while (NormOps.fastNormP2(f1) > rootTolerance) {
			CommonOps.invert(B, B_inv);
			CommonOps.mult(-1, B_inv, f1, dy);
			CommonOps.add(y1, dy, y2);

			rootFunctionEuler(t, xV, y2, step, f2);
			CommonOps.sub(f2, f1, df);

			CommonOps.mult(B, dy, tmpV);
			CommonOps.sub(df, tmpV, tmpV);
			CommonOps.multTransB(1 / SpecializedOps.elementSumSq(dy), tmpV, dy, tmpM);
			CommonOps.add(B, tmpM, B);

			DenseMatrix64F yTmp = y1;
			y1 = y2;
			y2 = yTmp;
			DenseMatrix64F fTmp = f1;
			f1 = f2;
			f2 = fTmp;
		}
		System.arraycopy(y1.getData(), 0, out, 0, ode.getDimensionOfVectorField());
	}

//	protected void rootFunction(double t, double[] x, double[] y, double step, double[] result) {
//		ode.computeVectorField(t, y, xDotTmp);
//		for (int i=0; i < ode.getDimensionOfVectorField(); i++)
//			result[i] = x[i] + step * xDotTmp[i] - y[i];
//	}
//
//	protected void findNextState(double t, double[] x, double step) {
//		Array2DRowRealMatrix Br = new Array2DRowRealMatrix(x.length, x.length);
//		for (int i=0; i < Br.getRowDimension(); i++)
//			for (int j=0; j < Br.getColumnDimension(); j++)
//				Br.setEntry(i, j, (i == j) ? 1.0 : 0.0);
//		RealMatrix B = Br;
//		ArrayRealVector y1 = new ArrayRealVector(x);
////		ArrayRealVector y2 = new ArrayRealVector(x2);
//		ArrayRealVector f1 = new ArrayRealVector(ode.getDimensionOfVectorField());
//		ArrayRealVector f2 = new ArrayRealVector(ode.getDimensionOfVectorField());
//		rootFunction(t, x, y1.getDataRef(), step, f1.getDataRef());
//		//		rootFunction(t, x, y2.getDataRef(), step, f2.getDataRef());
////		dy = y2.subtract(y1);
////		df = f2.subtract(f1);
//		do {
//
//			SingularValueDecomposition svd = new SingularValueDecomposition(B);
//			RealMatrix B_inv = svd.getSolver().getInverse();
//			RealVector dy = B_inv.operate(f1);
//			dy.mapMultiplyToSelf(-1);
//			ArrayRealVector y2 = y1.add(dy);
//
//			rootFunction(t, x, y2.getDataRef(), step, f2.getDataRef());
//			ArrayRealVector df = f2.subtract(f1);
//
//			RealMatrix Q = df.subtract(B.operate(dy)).outerProduct(dy);
//			double dysquare = dy.dotProduct(dy);
//			double dysquare_inv = 1 / dysquare;
//			for (int i=0; i < B.getRowDimension(); i++)
//				for (int j=0; j < B.getColumnDimension(); j++)
//					Q.multiplyEntry(i, j, dysquare_inv);
//			B = B.add(Q);
//
//			ArrayRealVector fTmp = f2;
//			f2 = f1;
//			f1 = fTmp;
//			ArrayRealVector vTmp = y2;
//			y2 = y1;
//			y1 = vTmp;
//			System.out.println("y1:" + y1.getEntry(0));
//		} while (f1.getNorm() > rootTolerance);
//		for (int i=0; i < y1.getDimension(); i++)
//			x[i] = y1.getEntry(i);
//	}

    public static void main(String args[]) {
    	MyOde ode = new MyOde();
    	MyEf ef = new MyEf();
    	MyObserver observer = new MyObserver();
    	MyEventObserver eventObserver = new MyEventObserver();
    	double step = 0.01;
    	ImplicitEulerSolver q = new ImplicitEulerSolver(step, 1e-4, 1e-4, 1e-4);
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
