package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

public class ArrayFiniteTrajectory extends AbstractFiniteTrajectory {

	protected double[] tSeries;
	protected double[][] xSeries;

	public ArrayFiniteTrajectory(double[] tSeries) {
		this.tSeries = tSeries;
		xSeries = new double[0][0];
	}

	public ArrayFiniteTrajectory(double[] tSeries, double[][] xSeries) {
		this.tSeries = tSeries;
		this.xSeries = xSeries;
	}

	public ArrayFiniteTrajectory(FiniteTrajectory tr) {
		tSeries = tr.gettSeries();
		xSeries = tr.getxSeries();
	}

	@Override
	public int getNumberOfStates() {
		return xSeries.length;
	}

	@Override
	public int getNumberOfTimePoints() {
		return tSeries.length;
	}

	@Override
	public double[] getxState(int timePointIndex) {
		double[] x = new double[getNumberOfStates()];
		for (int s=0; s < x.length; s++)
			x[s] = xSeries[s][timePointIndex];
		return x;
	}

	@Override
	public double getInitialtime() {
		return tSeries[0];
	}

	@Override
	public double getFinalTime() {
		return tSeries[tSeries.length - 1];
	}

	@Override
	public double[] gettSeries() {
		return tSeries.clone();
	}

	@Override
	public double[][] getxSeries() {
		return xSeries.clone();
	}

	@Override
	public double[] getxSeries(int s) {
		return xSeries[s].clone();
	}

	@Override
	public RealVector gettVector() {
		return new ArrayRealVector(gettSeries());
	}

	@Override
	public RealVector getxVector(int s) {
		return new ArrayRealVector(getxSeries(s));
	}

	@Override
	public List<RealVector> getxVectors() {
		List<RealVector> l = new ArrayList<RealVector>(getNumberOfStates());
		for (int s=0; s < getNumberOfStates(); s++)
			l.add(getxVector(s));
		return l;
	}

	@Override
	public Iterator<RealVector> xVectorIterator() {
		return getxVectors().iterator();
	}

	@Override
	protected int binarySearchTimePoints(double t) {
		return Arrays.binarySearch(tSeries, t);
	}

	@Override
	public double[] getLinearCombination(double[] coefficients) {
		return getLinearCombination(xSeries, coefficients);
	}

}
