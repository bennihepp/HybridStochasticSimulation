package ch.ethz.khammash.hybridstochasticsimulation.simulators.lsodar;

import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;
import ch.ethz.khammash.nativeode.TimepointIterator;

public class LsodarTimepointIterator<T extends ReactionNetworkModel> implements TimepointIterator {

	private double currentTime;
	private double[] tSeries;
	private int index;

	public LsodarTimepointIterator(double t0, double t1) {
		tSeries = new double[2];
		tSeries[0] = t0;
		tSeries[1] = t1;
		currentTime = tSeries[0];
		index = 0;
	}

	public LsodarTimepointIterator(double[] tSeries) {
		this.tSeries = tSeries;
		index = 0;
	}

	public LsodarTimepointIterator(FiniteTrajectoryRecorder<T> tr) {
		this.tSeries = tr.gettSeries();
	}

	@Override
	public double getInitialTimepoint() {
		return tSeries[0];
	}

	@Override
	public double getLastTimepoint() {
		return tSeries[tSeries.length - 1];
	}

	@Override
	public void reset() {
		index = 0;
	}

	@Override
	public double getCurrentTimepoint() {
		return currentTime;
	}

	public void setCurrentTime(double t) {
		this.currentTime = t;
		while ((index + 1) < tSeries.length && t > tSeries[index + 1])
			index++;
	}

	@Override
	public boolean hasNextTimepoint() {
		return (index + 1) < tSeries.length;
	}

	@Override
	public double getNextTimepoint() {
		while (currentTime >= tSeries[index])
			index++;
		currentTime = tSeries[index];
		return tSeries[index];
	}

}
