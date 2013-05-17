package ch.ethz.khammash.hybridstochasticsimulation;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class TrajectoryDistributionPlotData extends TrajectoryDistributionData implements PlotData {

	private DefaultPlotData plotData;

	public TrajectoryDistributionPlotData(RealVector tVector) {
		super(tVector);
		plotData = new DefaultPlotData(0);
	}

	public TrajectoryDistributionPlotData(RealVector tVector, RealMatrix xMeanMatrix, RealMatrix xStdDevMatrix) {
		super(tVector, xMeanMatrix, xStdDevMatrix);
		plotData = new DefaultPlotData(getNumberOfStates());
	}

	public TrajectoryDistributionPlotData(String[] names, double[] plotScales, RealVector tVector, RealMatrix xMeanMatrix,
			RealMatrix xStdDevMatrix) {
		super(tVector, xMeanMatrix, xStdDevMatrix);
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
	public void addState(RealVector xMeanVector, RealVector xStdDevVector) {
		addState(plotData.DEFAULT_NAME, 1.0, xMeanVector, xStdDevVector);
	}

	public void addState(String name, double plotScale, RealVector xMeanVector, RealVector xStdDevVector) {
		super.addState(xMeanVector, xStdDevVector);
		plotData.addState(name, plotScale);
	}

	@Override
	public void removeState(int s) {
		super.removeState(s);
		plotData.removeState(s);
	}

	public TrajectoryDistributionPlotData getSubsetData(int[] states) {
		return getSubsetData(states, null);
	}

	public TrajectoryDistributionPlotData getSubsetData(int[] states, double[] plotScales) {
		TrajectoryDistributionPlotData tdd = new TrajectoryDistributionPlotData(gettVector());
		for (int i=0; i < states.length; i++) {
			int s = states[i];
			double plotScale = (plotScales != null) ? plotScales[i] : getPlotScale(s);
			tdd.addState(getName(s), plotScale, getxMeanVector(s), getxStdDevVector(s));
		}
		return tdd;
	}

}
