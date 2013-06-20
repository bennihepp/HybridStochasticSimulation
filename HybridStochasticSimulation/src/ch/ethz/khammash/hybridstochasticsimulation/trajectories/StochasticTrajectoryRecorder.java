package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;


public class StochasticTrajectoryRecorder<T extends ReactionNetworkModel> implements TrajectoryRecorder<T> {

	private List<ReactionEvent> reactions;

	public StochasticTrajectoryRecorder() {
		reactions = new ArrayList<ReactionEvent>();
	}

	public int getNumberOfReactionEvents() {
		return reactions.size();
	}

	public ReactionEvent getReactionEvent(int i) {
		return reactions.get(i);
	}

	public Iterator<ReactionEvent> iterator() {
		return reactions.iterator();
	}

	@Override
	public void setModel(T model) {
	}

	@Override
	public void setInitialState(double t0, double[] x0) {
		setInitialState(t0, x0, x0.length);
	}

	@Override
	public void setInitialState(double t0, double[] x0, int numOfStates) {
		reportState(t0, x0);
//		handleReactionEvent(-1, t0, x0);
	}

	@Override
	public void setFinalState(double t1, double[] x1) {
		reportState(t1, x1);
//		handleReactionEvent(-1, t1, x1);
	}

	@Override
	public void reportState(double t, double[] x) {
		reactions.add(new ReactionEvent(t, x.clone()));
	}

//	@Override
//	public void handleReactionEvent(int reaction, double t, double[] newX) {
//		reactions.add(new ReactionEvent(reaction, t, newX.clone()));
//	}

	@Override
	public double[] getInterpolatedState(double t) {
		int index = findPreviousReactionEventIndex(t);
		ReactionEvent re = reactions.get(index);
		return re.getNewX().clone();
	}

	@Override
	public RealVector getInterpolatedStateVector(double t) {
		int index = findPreviousReactionEventIndex(t);
		ReactionEvent re = reactions.get(index);
		return new ArrayRealVector(re.getNewX());
	}

	@Override
	public double getInitialtime() {
		return reactions.get(0).getTime();
	}

	@Override
	public double getFinalTime() {
		return reactions.get(reactions.size()-1).getTime();
	}

	private int findPreviousReactionEventIndex(double t) {
		int index = Collections.binarySearch(reactions, t);
		if (index >= 0)
			return index;
		// See Java API
		int insertionPoint = -index - 1;
		if (insertionPoint == 0)
			throw new InvalidTimePointException("Time must be between start and end of simulation");
		return insertionPoint - 1;
	}

}
