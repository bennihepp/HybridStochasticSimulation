package ch.ethz.khammash.hybridstochasticsimulation;

import static com.google.common.base.Preconditions.*;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;


public class TrajectoryData {

	private RealVector tVector;
	private List<RealVector> xVectors;

	public TrajectoryData(RealVector tVector) {
		this.tVector = tVector;
		xVectors = new ArrayList<RealVector>();
	}

	public TrajectoryData(RealVector tVector, RealMatrix xMatrix) {
		checkArgument(tVector.getDimension() == xMatrix.getColumnDimension(),
				"Size of tVector must be equal to number of columns of xMatrix");
		this.tVector = tVector;
		xVectors = new ArrayList<RealVector>(xMatrix.getRowDimension());
		for (int s=0; s < xMatrix.getRowDimension(); s++)
			xVectors.add(xMatrix.getRowVector(s));
	}

	public TrajectoryData(RealVector tVector, List<RealVector> xVectors) {
		for (RealVector xVector : xVectors)
			checkArgument(tVector.getDimension() == xVector.getDimension(),
					"Size of tVector must be equal to size of each xVector");
		this.tVector = tVector;
		this.xVectors = new ArrayList<RealVector>(xVectors);
	}

	public TrajectoryData(TrajectoryData td) {
		this(td.gettVector(), td.getxVectors());
	}

	public RealVector gettVector() {
		return tVector;
	}

	public void settVector(RealVector tVector) {
		checkArgument(gettVector().getDimension() == tVector.getDimension(),
				"Size of tVector can not be changed");
		this.tVector = tVector;
	}

	public RealVector getxVector(int s) {
		checkElementIndex(s, getNumberOfStates(), "Expected 0<=s<getNumberOfStates()");
		return xVectors.get(s);
	}

	public void setxVector(int s, RealVector xVector) {
		checkElementIndex(s, getNumberOfStates(), "Expected 0<=s<getNumberOfStates()");
		xVectors.set(s, xVector);
	}

	public int getNumberOfStates() {
		return xVectors.size();
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

	public List<RealVector> getxVectors() {
		return xVectors;
	}

	public Iterator<RealVector> iterator() {
		return xVectors.iterator();
	}

	public SingleTrajectoryData getSingleTrajectoryData(int s) {
		checkElementIndex(s, getNumberOfStates(), "Expected 0<=s<getNumberOfStates()");
		return new SingleTrajectoryData(tVector, xVectors.get(s));
	}

	public TrajectoryData getSubsetData(int[] states) {
		TrajectoryData td = new TrajectoryData(gettVector());
		for (int s : states) {
			checkElementIndex(s, getNumberOfStates(), "Expected 0<=s<getNumberOfStates()");
			td.addState(getxVector(s));
		}
		return td;
	}

	public RealVector getLinearCombinationOfxVectors(double[] coefficients) {
		return getLinearCombinationOfVectors(xVectors, coefficients);
	}

	protected RealVector getLinearCombinationOfVectors(List<RealVector> vectors, double[] coefficients) {
		checkArgument(vectors.size() == coefficients.length, "Expected vectors.size() = coefficients.length but %s != %s",
				vectors.size(), coefficients.length);
		// TODO: check length of coefficients
		RealVector lc = new ArrayRealVector(vectors.get(0).getDimension());
		for (int s=0; s < vectors.size(); s++) {
			double coeff = coefficients[s];
			RealVector v = vectors.get(s);
			v = v.mapMultiply(coeff);
			lc = lc.add(v);
		}
		return lc;
	}

}
