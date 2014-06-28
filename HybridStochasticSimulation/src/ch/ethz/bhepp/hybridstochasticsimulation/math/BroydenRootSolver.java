package ch.ethz.bhepp.hybridstochasticsimulation.math;

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import org.ejml.ops.NormOps;
import org.ejml.ops.SpecializedOps;

public class BroydenRootSolver {

	protected MultivariateFunction f;
	protected double rootTolerance;
	protected int dimension;
	protected DenseMatrix64F B;
	protected DenseMatrix64F B_inv;
	protected DenseMatrix64F tmpM;
	protected DenseMatrix64F xV;
	protected DenseMatrix64F y1V;
	protected DenseMatrix64F y2V;
	protected DenseMatrix64F dyV;
	protected DenseMatrix64F q1V;
	protected DenseMatrix64F q2V;
	protected DenseMatrix64F dqV;
	protected DenseMatrix64F tmpV;

	public BroydenRootSolver(MultivariateFunction f) {
		this(f, 1e-6);
	}

	public BroydenRootSolver(MultivariateFunction f, double rootTolerance) {
		this.f = f;
		this.rootTolerance = rootTolerance;
		dimension = f.getDimension();
		B = new DenseMatrix64F(dimension, dimension);
		B_inv = new DenseMatrix64F(dimension, dimension);
		tmpM = new DenseMatrix64F(dimension, dimension);
		xV = new DenseMatrix64F(dimension, 1);
		y1V = new DenseMatrix64F(dimension, 1);
		y2V = new DenseMatrix64F(dimension, 1);
		dyV = new DenseMatrix64F(dimension, 1);
		q1V = new DenseMatrix64F(dimension, 1);
		q2V = new DenseMatrix64F(dimension, 1);
		dqV = new DenseMatrix64F(dimension, 1);
		tmpV = new DenseMatrix64F(dimension, 1);
	}

	// FIXME: Add maximum number of iterations to prevent getting stuck
	public double[] findRoot(double[] startingPoint) {
		double[] x = startingPoint;
		xV.setData(x);
		y1V.set(xV);

		double[] y1 = y1V.getData();
		double[] q1 = q1V.getData();
		f.computeValue(y1, q1);

		CommonOps.setIdentity(B);
		while (NormOps.fastNormP2(q1V) > rootTolerance) {
			CommonOps.invert(B, B_inv);
			CommonOps.mult(-1, B_inv, q1V, dyV);
			CommonOps.add(y1V, dyV, y2V);

			double[] y2 = y2V.getData();
			double[] q2 = q2V.getData();
			f.computeValue(y2, q2);
			CommonOps.sub(q2V, q1V, dqV);

			CommonOps.mult(B, dyV, tmpV);
			CommonOps.sub(dqV, tmpV, tmpV);
			double inverseElementSumOfDyV = 1 / SpecializedOps.elementSumSq(dyV);
			CommonOps.multTransB(inverseElementSumOfDyV, tmpV, dyV, tmpM);
			CommonOps.add(B, tmpM, B);

			DenseMatrix64F yTmp = y1V;
			y1V = y2V;
			y2V = yTmp;
			DenseMatrix64F fTmp = q1V;
			q1V = q2V;
			q2V = fTmp;
		}

		return y1V.getData().clone();
	}

}
