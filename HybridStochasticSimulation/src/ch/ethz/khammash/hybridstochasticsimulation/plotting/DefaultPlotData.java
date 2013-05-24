package ch.ethz.khammash.hybridstochasticsimulation.plotting;

import static com.google.common.base.Preconditions.*;

import java.util.ArrayList;
import java.util.List;


public class DefaultPlotData implements PlotData {

	protected final String DEFAULT_NAME = "<unnamed>";

	private List<String> names;
	private List<Double> plotScales;

	public DefaultPlotData(int numberOfStates) {
		names = new ArrayList<String>(numberOfStates);
		plotScales = new ArrayList<Double>(numberOfStates);
		for (int i=0; i < numberOfStates; i++) {
			names.add(DEFAULT_NAME);
			plotScales.add(Double.valueOf(1.0));
		}
	}

	public void removeState(int s) {
		checkElementIndex(s, getNumberOfStates(), "Expected 0<=s<getNumberOfStates()");
		names.remove(s);
		plotScales.remove(s);
	}

	public void addState(String name, double plotScale) {
		names.add(name);
		plotScales.add(Double.valueOf(plotScale));
	}

	@Override
	public String getName(int s) {
		checkElementIndex(s, getNumberOfStates(), "Expected 0<=s<getNumberOfStates()");
		return names.get(s);
	}

	public void setName(int s, String name) {
		checkElementIndex(s, getNumberOfStates(), "Expected 0<=s<getNumberOfStates()");
		names.set(s,  name);
	}

	@Override
	public String[] getNames() {
		return names.toArray(new String[0]);
	}

	public void setNames(String[] names) {
		checkArgument(names.length == getNumberOfStates(), "Expected names.length == getNumberOfStates()");
		this.names.clear();
		for (String name : names)
			this.names.add(name);
	}

	@Override
	public double getPlotScale(int s) {
		checkElementIndex(s, getNumberOfStates(), "Expected 0<=s<getNumberOfStates()");
		return plotScales.get(s).doubleValue();
	}

	public void setPlotScale(int s, double plotScale) {
		plotScales.set(s, Double.valueOf(plotScale));
	}

	@Override
	public double[] getPlotScales() {
		double[] d = new double[plotScales.size()];
		for (int i=0; i < plotScales.size(); i++)
			d[i] = plotScales.get(i);
		return d;
	}

	public void setPlotScales(double[] plotScales) {
		checkArgument(plotScales.length == getNumberOfStates(), "Expected names.length == getNumberOfStates()");
		this.plotScales.clear();
		for (double plotScale : plotScales)
			this.plotScales.add(plotScale);
	}

	@Override
	public int getNumberOfStates() {
		return names.size();
	}

	@Override
	public boolean isContinuous() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean isDiscrete() {
		// TODO Auto-generated method stub
		return false;
	}

}
