package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import java.util.Iterator;
import java.util.List;


public interface PlotData extends Trajectory {

	String getDescription();

	int getNumberOfStates();

	String getStateName(int state);

	List<String> getStateNames();

	Iterator<String> stateNameIterator();

	List<Double> getPlotScales();

	double getPlotScale(int state);

	boolean isContinuous();

	boolean isDiscrete();

	PlotData getSubsetData(Integer... states);

	PlotData getSubsetData(int[] states);

	PlotData getSubsetData(List<Integer> states);

	PlotData getSubsetData(int[] states, double[] plotScales);

	PlotData getSubsetData(List<Integer> states, List<Double> plotScales);

}
