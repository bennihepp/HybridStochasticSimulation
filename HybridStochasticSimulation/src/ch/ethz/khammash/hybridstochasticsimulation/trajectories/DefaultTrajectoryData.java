package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import static com.google.common.base.Preconditions.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;



public class DefaultTrajectoryData implements TrajectoryData {

	private RealVector tVector;
	private List<RealVector> xVectors;
	private List<String> stateNames;

	public DefaultTrajectoryData(RealVector tVector) {
		this.tVector = tVector;
		xVectors = new ArrayList<RealVector>();
		initStateNames();
	}

	public DefaultTrajectoryData(RealVector tVector, RealMatrix xMatrix) {
		checkArgument(tVector.getDimension() == xMatrix.getColumnDimension(),
				"Size of tVector must be equal to number of columns of xMatrix");
		this.tVector = tVector;
		xVectors = new ArrayList<RealVector>(xMatrix.getRowDimension());
		for (int s=0; s < xMatrix.getRowDimension(); s++)
			xVectors.add(xMatrix.getRowVector(s));
		initStateNames();
	}

	public DefaultTrajectoryData(RealVector tVector, List<RealVector> xVectors) {
		for (RealVector xVector : xVectors)
			checkArgument(tVector.getDimension() == xVector.getDimension(),
					"Size of tVector must be equal to size of each xVector");
		this.tVector = tVector;
		this.xVectors = new ArrayList<RealVector>(xVectors);
		initStateNames();
	}

	private void initStateNames() {
		stateNames = new ArrayList<String>(xVectors.size());
		for (int s=0; s < xVectors.size(); s++)
			stateNames.add("S" + s);
	}

	public DefaultTrajectoryData(DefaultTrajectoryData td) {
		this(td.gettVector(), td.getxVectors());
	}

	@Override
	public RealVector gettVector() {
		return tVector;
	}

	public void settVector(RealVector tVector) {
		checkArgument(gettVector().getDimension() == tVector.getDimension(),
				"Size of tVector can not be changed");
		this.tVector = tVector;
	}

	@Override
	public String getStateName(int s) {
		checkElementIndex(s, getNumberOfStates(), "Expected 0<=s<getNumberOfStates()");
		return stateNames.get(s);
	}

	public void setStateName(int s, String name) {
		checkElementIndex(s, getNumberOfStates(), "Expected 0<=s<getNumberOfStates()");
		stateNames.set(s, name);
	}

	@Override
	public RealVector getxVector(int s) {
		checkElementIndex(s, getNumberOfStates(), "Expected 0<=s<getNumberOfStates()");
		return xVectors.get(s);
	}

	public void setxVector(int s, RealVector xVector) {
		checkElementIndex(s, getNumberOfStates(), "Expected 0<=s<getNumberOfStates()");
		xVectors.set(s, xVector);
	}

	@Override
	public int getNumberOfStates() {
		return xVectors.size();
	}

	public void addState(String name, RealVector xVector) {
		checkArgument(gettVector().getDimension() == xVector.getDimension(),
				"Size of tVector must be equal to size of xVector");
		xVectors.add(xVector);
		stateNames.add(name);
	}

	public void addState(RealVector xVector) {
		addState("S" + stateNames.size(), xVector);
	}

	public void removeState(int s) {
		checkElementIndex(s, getNumberOfStates(), "Expected 0<=s<getNumberOfStates()");
		xVectors.remove(s);
		stateNames.remove(s);
	}

	@Override
	public List<RealVector> getxVectors() {
		return Collections.unmodifiableList(xVectors);
	}

	@Override
	public List<String> getStateNames() {
		return Collections.unmodifiableList(stateNames);
	}

	@Override
	public Iterator<RealVector> xVectorIterator() {
		return Collections.unmodifiableList(xVectors).iterator();
	}

	@Override
	public Iterator<String> stateNameIterator() {
		return Collections.unmodifiableList(stateNames).iterator();
	}

//	@Override
//	public DefaultSingleTrajectoryData getSingleTrajectoryData(int s) {
//		checkElementIndex(s, getNumberOfStates(), "Expected 0<=s<getNumberOfStates()");
//		return new DefaultSingleTrajectoryData(tVector, xVectors.get(s));
//	}

	@Override
	public DefaultTrajectoryData getSubsetData(int[] states) {
		DefaultTrajectoryData td = new DefaultTrajectoryData(gettVector());
		for (int s : states) {
			checkElementIndex(s, getNumberOfStates(), "Expected 0<=s<getNumberOfStates()");
			td.addState(getStateName(s), getxVector(s));
		}
		return td;
	}

	@Override
	public RealVector getLinearCombinationOfStates(double[] coefficients) {
		return getLinearCombinationOfStates(xVectors, coefficients);
	}

	@Override
	public RealVector getLinearCombinationOfStates(RealVector coefficients) {
		return getLinearCombinationOfStates(xVectors, coefficients.toArray());
	}

	protected RealVector getLinearCombinationOfStates(List<RealVector> vectors, double[] coefficients) {
		checkArgument(vectors.size() == coefficients.length, "Expected vectors.size() = coefficients.length but %s != %s",
				vectors.size(), coefficients.length);
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
