package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkElementIndex;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteDistributionPlotData;


public class VectorFiniteDistributionPlotData extends VectorFinitePlotData implements FiniteDistributionPlotData {

	private List<RealVector> xStdDevVectors;

	public VectorFiniteDistributionPlotData(double[] tSeries) {
		super(tSeries);
		initxStdDevVectors();
	}

	public VectorFiniteDistributionPlotData(double[] tSeries, double[][] xSeries, double[][] xStdDevSeries) {
		super(tSeries, xSeries);
		initxStdDevVectors(xStdDevSeries);
	}

	public VectorFiniteDistributionPlotData(RealVector tVector) {
		super(tVector);
		initxStdDevVectors();
	}

	public VectorFiniteDistributionPlotData(RealVector tVector, List<RealVector> xVectors, List<RealVector> xStdDevVectors) {
		super(tVector, xVectors);
		initxStdDevVectors(xStdDevVectors);
	}

	public VectorFiniteDistributionPlotData(RealVector tVector, RealMatrix xMatrix, RealMatrix xStdDevMatrix) {
		super(tVector, xMatrix);
		initxStdDevVectors(xStdDevMatrix);
	}

	public VectorFiniteDistributionPlotData(FiniteDistributionTrajectory tr) {
		super(tr);
		initxStdDevVectors(xStdDevVectors);
	}

	public VectorFiniteDistributionPlotData(FiniteDistributionPlotData td) {
		super(td);
		initxStdDevVectors(xStdDevVectors);
	}

	final private void initxStdDevVectors() {
		xStdDevVectors = new ArrayList<RealVector>();
	}

	final private void initxStdDevVectors(double[][] xStdDevSeries) {
		checkArgument(getNumberOfStates() == xStdDevSeries.length);
		for (int s=0; s < xStdDevSeries.length; s++)
			checkArgument(getNumberOfTimePoints() == xStdDevSeries[s].length);
		xStdDevVectors = new ArrayList<RealVector>(getNumberOfStates());
		for (int s=0; s < getNumberOfStates(); s++)
			xStdDevVectors.add(new ArrayRealVector(xStdDevSeries[s]));
	}

	final private void initxStdDevVectors(List<RealVector> xStdDevVectors) {
		checkArgument(getNumberOfStates() == xStdDevVectors.size());
		for (RealVector xVector : xStdDevVectors)
			checkArgument(getNumberOfTimePoints() == xVector.getDimension());
		xStdDevVectors = new ArrayList<RealVector>(getNumberOfStates());
		for (int s=0; s < getNumberOfStates(); s++)
			xStdDevVectors.add(xStdDevVectors.get(s));
	}

	final private void initxStdDevVectors(RealMatrix xStdDevMatrix) {
		checkArgument(getNumberOfStates() == xStdDevMatrix.getRowDimension());
		checkArgument(getNumberOfTimePoints() == xStdDevMatrix.getColumnDimension());
		xStdDevVectors = new ArrayList<RealVector>(getNumberOfStates());
		for (int s=0; s < getNumberOfStates(); s++)
			xStdDevVectors.add(xStdDevMatrix.getRowVector(s));
	}

	@Override
	public double[] getxStdDevState(int timePointIndex) {
		double[] x = new double[getNumberOfStates()];
		for (int s=0; s < x.length; s++)
			x[s] = xStdDevVectors.get(s).getEntry(timePointIndex);
		return x;
	}

	@Override
	public double[] getxStdDevSeries(int s) {
		return xStdDevVectors.get(s).toArray();
	}

	@Override
	public double[][] getxStdDevSeries() {
		double[][] xSeries = new double[getNumberOfStates()][getNumberOfTimePoints()];
		for (int s=0; s < getNumberOfStates(); s++)
			xSeries[s] = getxStdDevSeries(s);
		return xSeries;
	}

	@Override
	public RealVector getxStdDevVector(int s) {
		checkElementIndex(s, getNumberOfStates(), "Expected 0<=s<getNumberOfStates()");
		return xStdDevVectors.get(s);
	}

	public void setxStdDevVector(int s, RealVector xStdDevVector) {
		checkElementIndex(s, getNumberOfStates(), "Expected 0<=s<getNumberOfStates()");
		xStdDevVectors.set(s, xStdDevVector);
	}

	@Override
	public List<RealVector> getxStdDevVectors() {
		return Collections.unmodifiableList(xStdDevVectors);
	}

	@Override
	public Iterator<RealVector> xStdDevVectorIterator() {
		return Collections.unmodifiableList(xStdDevVectors).iterator();
	}

	public void addState(String name, RealVector xVector, RealVector xStdDevVector, double plotScale) {
		super.addState(name, xVector, plotScale);
		xStdDevVectors.add(xStdDevVector);
	}

	public void addState(String name, RealVector xVector, RealVector xStdDevVector) {
		super.addState(name, xVector);
		xStdDevVectors.add(xStdDevVector);
	}

	public void addState(RealVector xVector, RealVector xStdDevVector) {
		super.addState(xVector);
		xStdDevVectors.add(xStdDevVector);
	}

	@Override
	public void removeState(int s) {
		super.removeState(s);
		xStdDevVectors.remove(s);
	}

	@Override
	public VectorFiniteDistributionPlotData getSubsetData(Integer... states) {
		return getSubsetData(Arrays.asList(states));
	}

	@Override
	public VectorFiniteDistributionPlotData getSubsetData(List<Integer> states) {
		return getSubsetData(states, getPlotScales());
	}

	@Override
	public VectorFiniteDistributionPlotData getSubsetData(int[] states) {
		return getSubsetData(Arrays.asList(ArrayUtils.toObject(states)));
	}

	@Override
	public VectorFiniteDistributionPlotData getSubsetData(int[] states, double[] plotScales) {
		return getSubsetData(Arrays.asList(ArrayUtils.toObject(states)), Arrays.asList(ArrayUtils.toObject(plotScales)));
	}

	@Override
	public VectorFiniteDistributionPlotData getSubsetData(List<Integer> states, List<Double> plotScales) {
		VectorFiniteDistributionPlotData pd = new VectorFiniteDistributionPlotData(gettVector());
		pd.setDescription(getDescription());
		pd.setDiscrete(isDiscrete());
		for (int s : states) {
			checkElementIndex(s, getNumberOfStates(), "Expected 0<=s<getNumberOfStates()");
			pd.addState(getStateName(s), getxVector(s), getxStdDevVector(s), plotScales.get(s));
		}
		return pd;
	}

	@Override
	public double[] getLinearCombinationStdDev(double[] coefficients) {
		return getLinearCombination(xStdDevVectors, new ArrayRealVector(coefficients)).toArray();
	}

	@Override
	public RealVector getLinearCombinationStdDev(RealVector coefficients) {
		return getLinearCombination(xStdDevVectors, coefficients);
	}

}
