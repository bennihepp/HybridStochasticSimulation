package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import ch.ethz.khammash.hybridstochasticsimulation.simulators.ReactionEvent;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.ReactionEventHandler;


public class StochasticModelTrajectory implements ModelTrajectory, ReactionEventHandler {

	protected List<ReactionEvent> reactions;

	public StochasticModelTrajectory() {
		reactions = new ArrayList<ReactionEvent>();
	}

	@Override
	public void setInitialState(double t, double[] x) {
		handleReactionEvent(-1, t, x);
	}

	@Override
	public void setFinalState(double t, double[] x) {
		handleReactionEvent(-1, t, x);
	}

	@Override
	public void handleReactionEvent(int reaction, double t, double[] newX) {
		reactions.add(new ReactionEvent(reaction, t, newX.clone()));
	}

	private int getPreviousReactionEventIndex(double t) {
		int index = Collections.binarySearch(reactions, t);
		if (index >= 0)
			return index;
		// See Java API
		int insertionPoint = -index - 1;
		if (insertionPoint == 0)
			throw new InvalidTimePointException("Time must be between start and end of simulation");
		return insertionPoint - 1;
	}

	@Override
	public double[] getInterpolatedState(double t) {
		int index = getPreviousReactionEventIndex(t);
		ReactionEvent re = reactions.get(index);
		return re.getNewX().clone();
	}

	@Override
	public RealVector getInterpolatedStateVector(double t) {
		int index = getPreviousReactionEventIndex(t);
		ReactionEvent re = reactions.get(index);
		return new ArrayRealVector(re.getNewX());
	}

	public double getInitialtime() {
		return reactions.get(0).getTime();
	}

	public double getFinalTime() {
		return reactions.get(reactions.size()-1).getTime();
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

}
