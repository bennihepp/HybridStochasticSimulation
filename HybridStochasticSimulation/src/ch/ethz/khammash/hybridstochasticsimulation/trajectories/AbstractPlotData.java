package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import static com.google.common.base.Preconditions.checkElementIndex;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;


public abstract class AbstractPlotData implements PlotData {

	private final String DEFAULT_TITLE = "Trajectory";

	private String description;
	private List<String> stateNames;
	private List<Double> plotScales;
	private boolean isDiscrete;

	public AbstractPlotData() {
		this(0);
	}

	public AbstractPlotData(int numOfStates) {
		stateNames = new ArrayList<String>(numOfStates);
		for (int s=0; s < numOfStates; s++)
			stateNames.add("S" + s);
		description = DEFAULT_TITLE;
		initPlotScales();
		isDiscrete = false;
	}

	public AbstractPlotData(PlotData td) {
		stateNames = td.getStateNames();
		description = DEFAULT_TITLE;
		initPlotScales();
		isDiscrete = false;
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

	@Override
	public int getNumberOfStates() {
		return stateNames.size();
	}

	protected void addState(String name) {
		stateNames.add(name);
	}

	protected void addState() {
		addState("S" + stateNames.size());
	}

	public void removeState(int s) {
		stateNames.remove(s);
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
