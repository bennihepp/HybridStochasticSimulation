package ch.ethz.khammash.hybridstochasticsimulation;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class TrajectoryPlotData extends TrajectoryData implements PlotData {

	private DefaultPlotData plotData;

	public TrajectoryPlotData(RealVector tVector) {
		super(tVector);
		plotData = new DefaultPlotData(0);
	}

	public TrajectoryPlotData(RealVector tVector, RealMatrix xMatrix) {
		super(tVector, xMatrix);
		plotData = new DefaultPlotData(getNumberOfStates());
	}

	public TrajectoryPlotData(String[] names, double[] plotScales, RealVector tVector, RealMatrix xMatrix) {
		super(tVector, xMatrix);
		plotData = new DefaultPlotData(getNumberOfStates());
		for (int s=0; s < getNumberOfStates(); s++) {
			plotData.setName(s, names[s]);
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

	@Override
	public void addState(RealVector xVector) {
		addState(plotData.DEFAULT_NAME, 1.0, xVector);
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
		TrajectoryPlotData tdd = new TrajectoryPlotData(gettVector());
		for (int i=0; i < states.length; i++) {
			int s = states[i];
			double plotScale = (plotScales != null) ? plotScales[i] : getPlotScale(s);
			tdd.addState(getName(s), plotScale, getxVector(s));
		}
		return tdd;
	}

}
