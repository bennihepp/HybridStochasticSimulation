package ch.ethz.bhepp.hybridstochasticsimulation.math;

import org.apache.commons.math3.util.FastMath;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import org.ejml.ops.NormOps;
import org.ejml.ops.SpecializedOps;

import ch.ethz.bhepp.hybridstochasticsimulation.math.BroydenRootSolver;

public class BoundedBroydenRootSolver extends BroydenRootSolver {

	private double[] lowerBounds;
	private double[] upperBounds;

	public BoundedBroydenRootSolver(MultivariateFunction f, double[] lowerBounds, double[] upperBounds) {
		super(f);
		this.lowerBounds = lowerBounds;
		this.upperBounds = upperBounds;
	}

	public BoundedBroydenRootSolver(MultivariateFunction f, double rootTolerance, double[] lowerBounds, double[] upperBounds) {
		super(f, rootTolerance);
		this.lowerBounds = lowerBounds;
		this.upperBounds = upperBounds;
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
			double factor = Double.POSITIVE_INFINITY;
			for (int i=0; i < y2V.numRows; i++) {
				if (y2V.get(i) < lowerBounds[i]) {
					factor = FastMath.min(factor, FastMath.abs((y1V.get(i) - lowerBounds[i]) / dyV.get(i)));
				}
				if (y2V.get(i) > upperBounds[i]) {
					factor = FastMath.min(factor, FastMath.abs((upperBounds[i] - y1V.get(i)) / dyV.get(i)));
				}
			}
			if (factor < Double.POSITIVE_INFINITY) {
				CommonOps.add(y1V, factor / 2, dyV, y2V);
			}

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
