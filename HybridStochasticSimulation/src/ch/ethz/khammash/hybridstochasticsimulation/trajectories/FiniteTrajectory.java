package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import java.util.Arrays;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;


public class FiniteTrajectory implements Trajectory {

	protected double[] tSeries;
	protected double[][] xSeries;
	protected int index;

	public FiniteTrajectory(double[] tSeries) {
		this.tSeries = tSeries;
	}

	public double[] gettSeries() {
		return tSeries;
	}

	public double[][] getxSeries() {
		return xSeries;
	}

	protected void initialize(double[] x0) {
		initialize(x0, x0.length);
	}

	protected void initialize(double[] x0, int numberOfStates) {
		xSeries = new double[numberOfStates][tSeries.length];
		index = 0;
	}

	protected void setState(int index, double[] x) {
		for (int i=0; i < xSeries.length; i++)
			xSeries[i][index] = x[i];
	}

	private int getPreviousTimePointIndex(double t) {
		int index = Arrays.binarySearch(tSeries, t);
		if (index >= 0)
			return index;
		// See Java API
		int insertionPoint = -index - 1;
		if (insertionPoint == 0)
			throw new InvalidTimePointException("Time must be between start and end of simulation");
		return insertionPoint - 1;
	}

	// TODO: This is a trivial interpolation (the interpolated state is the state that was recorded before that timepoint)
	@Override
	public double[] getInterpolatedState(double t) {
		int index = getPreviousTimePointIndex(t);
		double[] x = new double[xSeries.length];
		for (int s=0; s < x.length; s++)
			x[s] = xSeries[s][index];
		return x;
	}

	// TODO: This is a trivial interpolation (the interpolated state is the state that was recorded before that timepoint)
	@Override
	public RealVector getInterpolatedStateVector(double t) {
		return new ArrayRealVector(getInterpolatedState(t));
	}

}
