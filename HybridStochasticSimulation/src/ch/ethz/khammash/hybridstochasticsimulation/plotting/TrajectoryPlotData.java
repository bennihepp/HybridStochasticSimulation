package ch.ethz.khammash.hybridstochasticsimulation.plotting;

import static com.google.common.base.Preconditions.*;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import ch.ethz.khammash.hybridstochasticsimulation.models.TrajectoryData;

public class TrajectoryPlotData extends TrajectoryData implements PlotData {

	private DefaultPlotData plotData;
	private boolean isDiscrete;

	public TrajectoryPlotData(RealVector tVector) {
		super(tVector);
		init();
	}

	public TrajectoryPlotData(RealVector tVector, RealMatrix xMatrix) {
		super(tVector, xMatrix);
		init();
	}

	public TrajectoryPlotData(String[] names, double[] plotScales, RealVector tVector, RealMatrix xMatrix) {
		super(tVector, xMatrix);
		init(names, plotScales);
	}

	private void init() {
		init(null, null);
	}

	private void init(String[] names, double[] plotScales) {
		plotData = new DefaultPlotData(getNumberOfStates());
		isDiscrete = false;
		if (names != null)
			checkArgument(names.length == getNumberOfStates(), "Expected names.length == getNumberOfStates()");
		if (plotScales != null)
			checkArgument(plotScales.length == getNumberOfStates(), "Expected plotScales.length == getNumberOfStates()");
		for (int s=0; s < getNumberOfStates(); s++) {
			if (names != null)
				plotData.setName(s, names[s]);
			if (plotScales != null)
				plotData.setPlotScale(s, Double.valueOf(plotScales[s]));
		}
	}

	@Override
	public String getName(int s) {
		return plotData.getName(s);
	}

	public void setName(int s, String name) {
		plotData.setName(s,  name);
	}

	@Override
	public String[] getNames() {
		return plotData.getNames();
	}

	public void setNames(String[] names) {
		plotData.setNames(names);
	}

	@Override
	public double getPlotScale(int s) {
		return plotData.getPlotScale(s);
	}

	public void setPlotScale(int s, double plotScale) {
		plotData.setPlotScale(s, Double.valueOf(plotScale));
	}

	@Override
	public double[] getPlotScales() {
		return plotData.getPlotScales();
	}

	public void setPlotScales(double[] plotScales) {
		plotData.setPlotScales(plotScales);
	}

	public void addState(String name, double plotScale, RealVector xVector) {
		super.addState(xVector);
		plotData.addState(name, plotScale);
	}

	public void removeState(int s) {
		super.removeState(s);
		plotData.removeState(s);
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
			tdd.addState(getName(s), plotScale, getxVector(s));
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
		return plotData.getTitle();
	}

	public void setTitle(String title) {
		plotData.setTitle(title);
	}

}
