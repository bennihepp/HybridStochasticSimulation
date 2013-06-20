package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

public abstract class AbstractFiniteTrajectory implements FiniteTrajectory {

	protected int index;

	// TODO: This is a trivial interpolation (the interpolated state is the state that was recorded before that timepoint)
	@Override
	public double[] getInterpolatedState(double t) {
		int index = findPreviousTimePointIndex(t);
		return getxState(index);
	}

	// TODO: This is a trivial interpolation (the interpolated state is the state that was recorded before that timepoint)
	@Override
	public RealVector getInterpolatedStateVector(double t) {
		return new ArrayRealVector(getInterpolatedState(t));
	}

	protected int binarySearchTimePoints(double t) {
		return Arrays.binarySearch(gettSeries(), t);
	}

	private int findPreviousTimePointIndex(double t) {
		int index = binarySearchTimePoints(t);
		if (index >= 0)
			return index;
		// See Java API
		int insertionPoint = -index - 1;
		if (insertionPoint == 0)
			throw new InvalidTimePointException("Time must be between start and end of simulation");
		return insertionPoint - 1;
	}

	@Override
	public double[] getLinearCombination(double[] coefficients) {
		return getLinearCombination(getxSeries(), coefficients);
	}

	@Override
	public RealVector getLinearCombination(RealVector coefficients) {
		return getLinearCombination(getxVectors(), coefficients);
	}

	protected static double[] getLinearCombination(double[][] xSeries, double[] coefficients) {
		checkArgument(xSeries.length == coefficients.length, "Expected xSeries.length = coefficients.length but %s != %s",
				xSeries.length, coefficients.length);
		double[] lc = new double[xSeries[0].length];
		for (int s=0; s < xSeries.length; s++) {
			double coeff = coefficients[s];
			double[] x = xSeries[s];
			for (int i=0; i < lc.length; i++)
				lc[i] += coeff * x[i];
		}
		return lc;
	}

	protected static RealVector getLinearCombination(List<RealVector> vectors, RealVector coefficients) {
		checkArgument(vectors.size() == coefficients.getDimension(),
				"Expected vectors.size() = coefficients.getDimension but %s != %s",
				vectors.size(), coefficients.getDimension());
		RealVector lc = new ArrayRealVector(vectors.get(0).getDimension());
		for (int s=0; s < vectors.size(); s++) {
			double coeff = coefficients.getEntry(s);
			RealVector v = vectors.get(s);
			v = v.mapMultiply(coeff);
			lc = lc.add(v);
		}
		return lc;
	}

}
