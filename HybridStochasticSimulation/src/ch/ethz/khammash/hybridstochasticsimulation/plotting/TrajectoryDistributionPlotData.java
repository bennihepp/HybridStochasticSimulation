package ch.ethz.khammash.hybridstochasticsimulation.plotting;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkElementIndex;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;


public class TrajectoryDistributionPlotData extends TrajectoryPlotData {

	private List<RealVector> xStdDevVectors;

	public TrajectoryDistributionPlotData(RealVector tVector) {
		super(tVector);
		initxStdDevVectors();
	}

	public TrajectoryDistributionPlotData(RealVector tVector, RealMatrix xMeanMatrix, RealMatrix xStdDevMatrix) {
		super(tVector, xMeanMatrix);
		initxStdDevVectors(xStdDevMatrix);
	}

	public TrajectoryDistributionPlotData(String[] names, RealVector tVector, RealMatrix xMeanMatrix,
			RealMatrix xStdDevMatrix) {
		super(names,  tVector, xMeanMatrix);
		initxStdDevVectors(xStdDevMatrix);
	}

	public TrajectoryDistributionPlotData(String[] names, double[] plotScales, RealVector tVector, RealMatrix xMeanMatrix,
			RealMatrix xStdDevMatrix) {
		super(names, plotScales, tVector, xMeanMatrix);
		initxStdDevVectors(xStdDevMatrix);
	}

	private void initxStdDevVectors() {
		xStdDevVectors = new ArrayList<RealVector>();
	}

	private void initxStdDevVectors(RealMatrix xStdDevMatrix) {
		checkArgument(gettVector().getDimension() == xStdDevMatrix.getColumnDimension(),
				"Size of tVector must be equal to number of columns of xStdDevMatrix");
		xStdDevVectors = new ArrayList<RealVector>(getNumberOfStates());
		for (int s=0; s < xStdDevMatrix.getRowDimension(); s++)
			xStdDevVectors.add(xStdDevMatrix.getRowVector(s));
	}

	public RealVector getxMeanVector(int s) {
		return getxVector(s);
	}

	public void setxMeanVector(int s, RealVector xMeanVector) {
		setxVector(s, xMeanVector);
	}

	public RealVector getxStdDevVector(int s) {
		checkElementIndex(s, getNumberOfStates(), "Expected 0<=s<getNumberOfStates()");
		return xStdDevVectors.get(s);
	}

	public void setxStdDevVector(int s, RealVector xStdDevVector) {
		checkElementIndex(s, getNumberOfStates(), "Expected 0<=s<getNumberOfStates()");
		xStdDevVectors.set(s, xStdDevVector);
	}

	public List<RealVector> getxMeanVectors() {
		return getxVectors();
	}

	public List<RealVector> getxStdDevVectors() {
		return Collections.unmodifiableList(xStdDevVectors);
	}

	public Iterator<RealVector> xMeanIterator() {
		return xVectorIterator();
	}

	public Iterator<RealVector> xStdDevIterator() {
		return Collections.unmodifiableList(xStdDevVectors).iterator();
	}

	@Override
	public void addState(String name, double plotScale, RealVector xMeanVector) {
		RealVector xStdDevVector = new ArrayRealVector(xMeanVector.getDimension());
		addState(name, plotScale, xMeanVector, xStdDevVector);
	}

	public void addState(String name, double plotScale, RealVector xMeanVector, RealVector xStdDevVector) {
		super.addState(name, plotScale, xMeanVector);
		xStdDevVectors.add(xStdDevVector);
	}

	@Override
	public void removeState(int s) {
		super.removeState(s);
		xStdDevVectors.remove(s);
	}

	@Override
	public TrajectoryDistributionPlotData getSubsetData(int[] states) {
		return getSubsetData(states, null);
	}

	@Override
	public TrajectoryDistributionPlotData getSubsetData(int[] states, double[] plotScales) {
		if (plotScales != null)
			checkArgument(states.length == plotScales.length, "Expected states.length == plotScales.length");
		TrajectoryDistributionPlotData tdd = new TrajectoryDistributionPlotData(gettVector());
		tdd.setDiscrete(isDiscrete());
		for (int i=0; i < states.length; i++) {
			int s = states[i];
			checkElementIndex(s, getNumberOfStates(), "Expected 0<=s<getNumberOfStates()");
			double plotScale = (plotScales != null) ? plotScales[i] : getPlotScale(s);
			tdd.addState(getStateName(s), plotScale, getxMeanVector(s), getxStdDevVector(s));
		}
		return tdd;
	}

	public RealVector getLinearCombinationOfxMeanVectors(double[] coefficients) {
		return super.getLinearCombinationOfStates(coefficients);
	}

	public RealVector getLinearCombinationOfxStdDevVectors(double[] coefficients) {
		return getLinearCombinationOfStates(xStdDevVectors, coefficients);
	}

}
