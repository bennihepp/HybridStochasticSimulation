package ch.ethz.bhepp.hybridstochasticsimulation.trajectories;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkElementIndex;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class VectorFiniteTrajectory extends AbstractFiniteTrajectory {

	private static final long serialVersionUID = -55370111917885210L;

	private RealVector tVector;
	private List<RealVector> xVectors;

	public VectorFiniteTrajectory(double[] tSeries) {
		this.tVector = new ArrayRealVector(tSeries);
		xVectors = new ArrayList<RealVector>(0);
	}

	public VectorFiniteTrajectory(double[] tSeries, double[][] xSeries) {
		this.tVector = new ArrayRealVector(tSeries);
		xVectors = new ArrayList<RealVector>(xSeries.length);
		for (int s=0; s < xSeries.length; s++)
			xVectors.add(new ArrayRealVector(xSeries[s]));
	}

	public VectorFiniteTrajectory(RealVector tVector) {
		this.tVector = tVector.copy();
		xVectors = new ArrayList<RealVector>(0);
	}

	public VectorFiniteTrajectory(RealVector tVector, List<RealVector> xVectors) {
		for (RealVector xVector : xVectors)
			checkArgument(tVector.getDimension() == xVector.getDimension(),
					"Size of tVector must be equal to size of each xVector");
		this.tVector = tVector.copy();
		this.xVectors = copyVectorList(xVectors);
	}

	public VectorFiniteTrajectory(RealVector tVector, RealMatrix xMatrix) {
		checkArgument(tVector.getDimension() == xMatrix.getColumnDimension(),
				"Size of tVector must be equal to number of columns of xMatrix");
		this.tVector = tVector.copy();
		xVectors = new ArrayList<RealVector>(xMatrix.getRowDimension());
		for (int s=0; s < xMatrix.getRowDimension(); s++)
			xVectors.add(xMatrix.getRowVector(s));
	}

	public VectorFiniteTrajectory(FiniteTrajectory tr) {
		tVector = tr.gettVector().copy();
		xVectors = new ArrayList<RealVector>(tr.getNumberOfStates());
		for (int s=0; s < tr.getNumberOfStates(); s++)
			xVectors.add(tr.getxVector(s));
	}

	final private List<RealVector> copyVectorList(List<RealVector> src) {
		List<RealVector> dst = new ArrayList<RealVector>(src.size());
		for (RealVector vector : src)
			dst.add(new ArrayRealVector(vector));
		return dst;
	}

	@Override
	public int getNumberOfStates() {
		return xVectors.size();
	}

	@Override
	public int getNumberOfTimePoints() {
		return tVector.getDimension();
	}

	@Override
	public double[] getxState(int timePointIndex) {
		double[] x = new double[getNumberOfStates()];
		for (int s=0; s < x.length; s++)
			x[s] = xVectors.get(s).getEntry(timePointIndex);
		return x;
	}

	@Override
	public double getxState(int state, int timePointIndex) {
		return xVectors.get(state).getEntry(timePointIndex);
	}

	@Override
	public double getInitialtime() {
		return tVector.getEntry(0);
	}

	@Override
	public double getFinalTime() {
		return tVector.getEntry(tVector.getDimension() - 1);
	}

	@Override
	public double[] gettSeries() {
		return tVector.toArray();
	}

	@Override
	public double[] getxSeries(int s) {
		return xVectors.get(s).toArray();
	}

	@Override
	public double[][] getxSeries() {
		double[][] xSeries = new double[getNumberOfStates()][getNumberOfTimePoints()];
		for (int s=0; s < getNumberOfStates(); s++)
			xSeries[s] = getxSeries(s);
		return xSeries;
	}

	@Override
	public RealVector gettVector() {
		return tVector.copy();
	}

	@Override
	public RealVector getxVector(int s) {
		return xVectors.get(s).copy();
	}

	@Override
	public List<RealVector> getxVectors() {
		return Collections.unmodifiableList(xVectors);
	}

	@Override
	public Iterator<RealVector> xVectorIterator() {
		return Collections.unmodifiableList(xVectors).iterator();
	}

	public void addState(RealVector xVector) {
		checkArgument(gettVector().getDimension() == xVector.getDimension(),
				"Size of tVector must be equal to size of xVector");
		xVectors.add(xVector);
	}

	public void removeState(int s) {
		checkElementIndex(s, getNumberOfStates(), "Expected 0<=s<getNumberOfStates()");
		xVectors.remove(s);
	}

	@Override
	public RealVector getLinearCombination(RealVector coefficients) {
		return getLinearCombination(xVectors, coefficients);
	}

}
