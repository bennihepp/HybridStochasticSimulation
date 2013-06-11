package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import static com.google.common.base.Preconditions.checkArgument;
import static com.google.common.base.Preconditions.checkElementIndex;

import java.util.ArrayList;
import java.util.List;


import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;



public class TrajectoryPlotData extends DefaultTrajectoryData implements PlotData {

	private final String DEFAULT_TITLE = "Trajectory";

	private String title;
	private List<Double> plotScales;
	private boolean isDiscrete;

	public TrajectoryPlotData(RealVector tVector) {
		super(tVector);
		title = DEFAULT_TITLE;
		initPlotScales();
		isDiscrete = false;
	}

	public TrajectoryPlotData(RealVector tVector, RealMatrix xMatrix) {
		super(tVector, xMatrix);
		title = DEFAULT_TITLE;
		initPlotScales();
		isDiscrete = false;
	}

	public TrajectoryPlotData(String[] names, RealVector tVector, RealMatrix xMatrix) {
		super(tVector, xMatrix);
		title = DEFAULT_TITLE;
		initStateNames(names);
		initPlotScales();
		isDiscrete = false;
	}

	public TrajectoryPlotData(String[] names, double[] plotScales, RealVector tVector, RealMatrix xMatrix) {
		super(tVector, xMatrix);
		title = DEFAULT_TITLE;
		initStateNames(names);
		initPlotScales(plotScales);
		isDiscrete = false;
	}

	private void initStateNames(String[] names) {
		checkArgument(names.length == getNumberOfStates(), "Expected names.length == getNumberOfStates()");
		for (int s=0; s < getNumberOfStates(); s++)
			setStateName(s, names[s]);
	}

	private void initPlotScales() {
		plotScales = new ArrayList<Double>(getNumberOfStates());
		for (int s=0; s < getNumberOfStates(); s++)
			plotScales.add(s, 1.0);
	}

	private void initPlotScales(double[] plotScales) {
		checkArgument(plotScales.length == getNumberOfStates(), "Expected plotScales.length == getNumberOfStates()");
		this.plotScales = new ArrayList<Double>(plotScales.length);
		for (int s=0; s < getNumberOfStates(); s++)
			this.plotScales.add(s, plotScales[s]);
	}

	@Override
	public double getPlotScale(int s) {
		return plotScales.get(s);
	}

	public void setPlotScale(int s, double plotScale) {
		plotScales.set(s, plotScale);
	}

	public void setPlotScales(double[] plotScales) {
		checkArgument(plotScales.length == getNumberOfStates(), "Expected plotScales.length == getNumberOfStates()");
		for (int s=0; s < getNumberOfStates(); s++)
			setPlotScale(s, plotScales[s]);
	}

	public void addState(String name, double plotScale, RealVector xVector) {
		super.addState(name, xVector);
		plotScales.add(plotScale);
	}

	@Override
	public void addState(String name, RealVector xVector) {
		super.addState(name, xVector);
		plotScales.add(1.0);
	}

	@Override
	public void addState(RealVector xVector) {
		super.addState(xVector);
		plotScales.add(1.0);
	}

	@Override
	public void removeState(int s) {
		super.removeState(s);
		plotScales.remove(s);
	}

	public TrajectoryPlotData getSubsetData(int[] states) {
		return getSubsetData(states, null);
	}

	public TrajectoryPlotData getSubsetData(int[] states, double[] plotScales) {
		if (plotScales != null)
			checkArgument(states.length == plotScales.length, "Expected states.length == plotScales.length");
		TrajectoryPlotData tdd = new TrajectoryPlotData(gettVector());
		tdd.setDiscrete(isDiscrete());
		for (int i=0; i < states.length; i++) {
			int s = states[i];
			checkElementIndex(s, getNumberOfStates(), "Expected 0<=s<getNumberOfStates()");
			double plotScale = (plotScales != null) ? plotScales[i] : getPlotScale(s);
			tdd.addState(getStateName(s), plotScale, getxVector(s));
		}
		return tdd;
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

	@Override
	public String getTitle() {
		return title;
	}

	public void setTitle(String title) {
		this.title = title;
	}

}
