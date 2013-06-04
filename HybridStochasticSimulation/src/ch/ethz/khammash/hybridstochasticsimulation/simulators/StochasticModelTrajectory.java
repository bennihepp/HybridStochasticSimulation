package ch.ethz.khammash.hybridstochasticsimulation.simulators;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.math3.linear.RealVector;

import ch.ethz.khammash.hybridstochasticsimulation.models.ModelTrajectory;


public class StochasticModelTrajectory implements ModelTrajectory, ReactionEventHandler {

	protected List<ReactionEvent> reactions;

	public StochasticModelTrajectory() {
		reactions = new ArrayList<ReactionEvent>();
	}

	@Override
	public void setInitialState(double t0, double[] x0) {
		handleReactionEvent(-1, t0, x0);
	}

	@Override
	public void setFinalState(double t1, double[] x1) {
		handleReactionEvent(-1, t1, x1);
	}

	@Override
	public void handleReactionEvent(int reaction, double t, double[] newX) {
		reactions.add(new ReactionEvent(reaction, t, newX.clone()));
	}

	public double[] getInterpolatedState(double t) {
		int index = Collections.binarySearch(reactions, t);
		if (index >= 0)
			return reactions.get(index).getNewX();
		// See Java API
		int insertionPoint = -index - 1;
		if (insertionPoint == 0)
			throw new InvalidTimePointException("Time must be between start and end of simulation");
		ReactionEvent re = reactions.get(insertionPoint - 1);
		return re.getNewX();
	}

	public RealVector getInterpolatedStateVector(double t) {
		int index = Collections.binarySearch(reactions, t);
		if (index >= 0)
			return reactions.get(index).getNewXVector();
		// See Java API
		int insertionPoint = -index - 1;
		if (insertionPoint == 0)
			throw new InvalidTimePointException("Time must be between start and end of simulation");
		ReactionEvent re = reactions.get(insertionPoint - 1);
		return re.getNewXVector();
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

//	public int getIndexOfReaction(int i) {
//		return reactions.get(i).getReaction();
//	}
//
//	public double getTimeOfReaction(int i) {
//		return reactions.get(i).getTime();
//	}
//
//	public double[] getStateAfterReaction(int i) {
//		return reactions.get(i).getNewX();
//	}

	public ReactionEvent getReactionEvent(int i) {
		return reactions.get(i);
	}

	public Iterator<ReactionEvent> iterator() {
		return reactions.iterator();
	}

}
