package ch.ethz.khammash.hybridstochasticsimulation.simulators;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;


public class StochasticModelTrajectory implements ModelTrajectory, ReactionHandler {

	protected List<ReactionRecord> reactions;

	public StochasticModelTrajectory() {
		reactions = new ArrayList<ReactionRecord>();
	}

	@Override
	public void setInitialState(double t, double[] x) {
		handleReaction(-1, t, x);
	}

	@Override
	public void setFinalState(double t, double[] x) {
		handleReaction(-1, t, x);
	}

	@Override
	public void handleReaction(int reaction, double t, double[] newX) {
		reactions.add(new ReactionRecord(reaction, t, newX.clone()));
	}

	public double[] getInterpolatedState(double t) {
		int index = Collections.binarySearch(reactions, t);
		if (index >= 0)
			return reactions.get(index).getNewX();
		// See Java API
		int insertionPoint = -index - 1;
		if (insertionPoint == 0)
			throw new InvalidTimePointException("Time must be between start and end of simulation");
		ReactionRecord rr = reactions.get(insertionPoint - 1);
		return rr.getNewX();
	}

	public double getInitialtime() {
		return reactions.get(0).getTime();
	}

	public double getFinalTime() {
		return reactions.get(reactions.size()-1).getTime();
	}

	public int getNumberOfReactions() {
		return reactions.size();
	}

	public int getIndexOfReaction(int i) {
		return reactions.get(i).getReaction();
	}

	public double getTimeOfReaction(int i) {
		return reactions.get(i).getTime();
	}

	public double[] getStateAfterReaction(int i) {
		return reactions.get(i).getNewX();
	}

	public ReactionRecord getRecordOfReaction(int i) {
		return reactions.get(i);
	}

}
