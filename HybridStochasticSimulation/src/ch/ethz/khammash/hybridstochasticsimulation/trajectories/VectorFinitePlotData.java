package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkElementIndex;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;


public class VectorFinitePlotData extends VectorFiniteTrajectory implements FinitePlotData {

	private static final long serialVersionUID = -7723258694878731194L;

	private final String DEFAULT_TITLE = "Trajectory";

	private String description;
	private boolean isDiscrete;
	private List<String> stateNames;
	private List<Double> plotScales;

	public VectorFinitePlotData(double[] tSeries) {
		super(tSeries);
		description = DEFAULT_TITLE;
		isDiscrete = false;
		initStateNames();
		initPlotScales();
	}

	public VectorFinitePlotData(double[] tSeries, double[][] xSeries) {
		super(tSeries, xSeries);
		description = DEFAULT_TITLE;
		isDiscrete = false;
		initStateNames();
		initPlotScales();
	}

	public VectorFinitePlotData(RealVector tVector) {
		super(tVector);
		description = DEFAULT_TITLE;
		isDiscrete = false;
		initStateNames();
		initPlotScales();
	}

	public VectorFinitePlotData(RealVector tVector, List<RealVector> xVectors) {
		super(tVector, xVectors);
		description = DEFAULT_TITLE;
		isDiscrete = false;
		initStateNames();
		initPlotScales();
	}

	public VectorFinitePlotData(RealVector tVector, RealMatrix xMatrix) {
		super(tVector, xMatrix);
		description = DEFAULT_TITLE;
		isDiscrete = false;
		initStateNames();
		initPlotScales();
	}

	public VectorFinitePlotData(FiniteTrajectory tr) {
		super(tr);
		description = DEFAULT_TITLE;
		isDiscrete = false;
		initStateNames();
		initPlotScales();
	}

	public VectorFinitePlotData(FinitePlotData td) {
		super(td);
		description = td.getDescription();
		isDiscrete = td.isDiscrete();
		stateNames = new ArrayList<String>(td.getStateNames());
		plotScales = new ArrayList<Double>(td.getPlotScales());
		description = DEFAULT_TITLE;
		isDiscrete = false;
	}

	final private void initStateNames() {
		stateNames = new ArrayList<String>(getNumberOfStates());
		for (int s=0; s < getNumberOfStates(); s++)
			stateNames.add("S" + s);
	}

	final private void initPlotScales() {
		plotScales = new ArrayList<Double>(getNumberOfStates());
		for (int s=0; s < getNumberOfStates(); s++)
			plotScales.add(s, 1.0);
	}

	@Override
	public String getDescription() {
		return description;
	}

	public void setDescription(String description) {
		this.description = description;
	}

	@Override
	public double getPlotScale(int s) {
		return plotScales.get(s);
	}

	public void setPlotScale(int s, double plotScale) {
		plotScales.set(s, plotScale);
	}

	public void setPlotScales(Double... plotScales) {
		setPlotScales(Arrays.asList(plotScales));
	}

	public void setPlotScales(double[] plotScales) {
		checkArgument(plotScales.length == getNumberOfStates(), "Expected plotScales.length == getNumberOfStates()");
		for (int s=0; s < getNumberOfStates(); s++)
			setPlotScale(s, plotScales[s]);
	}

	public void setPlotScales(List<Double> plotScales) {
		checkArgument(plotScales.size() == getNumberOfStates(), "Expected plotScales.length == getNumberOfStates()");
		for (int s=0; s < getNumberOfStates(); s++)
			setPlotScale(s, plotScales.get(s));
	}

	@Override
	public List<Double> getPlotScales() {
		return Collections.unmodifiableList(plotScales);
	}

	@Override
	public List<String> getStateNames() {
		return Collections.unmodifiableList(stateNames);
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

	public void setStateNames(String[] names) {
		checkArgument(names.length == getNumberOfStates(), "Expected names.length == getNumberOfStates()");
		for (int s=0; s < getNumberOfStates(); s++)
			setStateName(s, names[s]);
	}

	public void setStateNames(List<String> names) {
		checkArgument(names.size() == getNumberOfStates(), "Expected names.length == getNumberOfStates()");
		for (int s=0; s < getNumberOfStates(); s++)
			setStateName(s, names.get(s));
	}

	public void addState(String name, RealVector xVector, double plotScale) {
		super.addState(xVector);
		stateNames.add(name);
		plotScales.add(plotScale);
	}

	public void addState(String name, RealVector xVector) {
		addState(name, xVector, 1.0);
	}

	public void addState(RealVector xVector) {
		addState("S" + getNumberOfStates(), xVector);
	}

	public void removeState(int s) {
		checkElementIndex(s, getNumberOfStates(), "Expected 0<=s<getNumberOfStates()");
		super.removeState(s);
		stateNames.remove(s);
		plotScales.remove(s);
	}

//	@Override
//	public DefaultSingleTrajectoryData getSingleTrajectoryData(int s) {
//		checkElementIndex(s, getNumberOfStates(), "Expected 0<=s<getNumberOfStates()");
//		return new DefaultSingleTrajectoryData(tVector, xVectors.get(s));
//	}

	@Override
	public VectorFinitePlotData getSubsetData(Integer... states) {
		return getSubsetData(Arrays.asList(states));
	}

	@Override
	public VectorFinitePlotData getSubsetData(List<Integer> states) {
		return getSubsetData(states, plotScales);
	}

	@Override
	public VectorFinitePlotData getSubsetData(int[] states) {
		return getSubsetData(Arrays.asList(ArrayUtils.toObject(states)));
	}

	@Override
	public VectorFinitePlotData getSubsetData(int[] states, double[] plotScales) {
		return getSubsetData(Arrays.asList(ArrayUtils.toObject(states)), Arrays.asList(ArrayUtils.toObject(plotScales)));
	}

	@Override
	public VectorFinitePlotData getSubsetData(List<Integer> states, List<Double> plotScales) {
		VectorFinitePlotData pd = new VectorFinitePlotData(gettVector());
		pd.setDescription(getDescription());
		pd.setDiscrete(isDiscrete());
		for (int i=0; i < states.size(); i++) {
			int s = states.get(i);
			double plotScale = plotScales.get(i);
			checkElementIndex(s, getNumberOfStates(), "Expected 0<=s<getNumberOfStates()");
			pd.addState(getStateName(s), getxVector(s), plotScale);
		}
		return pd;
	}

	@Override
	public Iterator<String> stateNameIterator() {
		return Collections.unmodifiableList(stateNames).iterator();
	}

	public void setContinuous() {
		isDiscrete = false;
	}

	public void setDiscrete() {
		isDiscrete = true;
	}

	public void setDiscrete(boolean isDiscrete) {
		this.isDiscrete = isDiscrete;
	}

	@Override
	public boolean isContinuous() {
		return !isDiscrete;
	}

	@Override
	public boolean isDiscrete() {
		return isDiscrete;
	}

}
