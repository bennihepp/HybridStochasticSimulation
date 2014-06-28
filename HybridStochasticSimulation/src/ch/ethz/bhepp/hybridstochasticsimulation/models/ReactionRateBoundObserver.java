package ch.ethz.bhepp.hybridstochasticsimulation.models;

import ch.ethz.bhepp.hybridstochasticsimulation.simulators.PDMPEventObserver;


public class ReactionRateBoundObserver implements PDMPEventObserver, OptionalEventObserver {

	public static enum BoundType {
		NONE, LOWER, UPPER, BOTH,
	}

	private ReactionRateBoundEventListener listener;
	private StochasticReactionNetworkModel model;
	private int reaction;
	private double lowerBound;
	private double upperBound;
	private BoundType boundType;

	public ReactionRateBoundObserver(ReactionRateBoundEventListener listener, StochasticReactionNetworkModel model, int reaction) {
		this.listener = listener;
		this.model = model;
		this.reaction = reaction;
		this.boundType = BoundType.NONE;
	}

	public ReactionRateBoundObserver(ReactionRateBoundEventListener listener, StochasticReactionNetworkModel model, int reaction, double bound, BoundType boundType) {
		this.listener = listener;
		this.model = model;
		this.reaction = reaction;
		lowerBound = bound;
		upperBound = bound;
		this.boundType = boundType;
	}

	public ReactionRateBoundObserver(ReactionRateBoundEventListener listener, StochasticReactionNetworkModel model, int reaction, double min, double max) {
		this.listener = listener;
		this.model = model;
		this.reaction = reaction;
		this.lowerBound = min;
		this.upperBound = max;
		this.boundType = BoundType.BOTH;
	}

	public double getLowerBound() {
		return lowerBound;
	}

	public void setLowerBound(double lowerBound) {
		this.lowerBound = lowerBound;
	}

	public double getUpperBound() {
		return upperBound;
	}

	public void setUpperBound(double upperBound) {
		this.upperBound = upperBound;
	}

	public BoundType getBoundType() {
		return boundType;
	}

	public void setBoundType(BoundType boundType) {
		this.boundType = boundType;
	}

	final public boolean isOutOfBounds(double t, double[] x) {
		double rate = model.computePropensity(reaction, t, x);
		switch (boundType) {
		case LOWER:
			return rate <= lowerBound;
		case UPPER:
			return rate >= upperBound;
		case BOTH:
			return (rate <= lowerBound) || (rate >= upperBound);
		case NONE:
		default:
			return false;
		}
	}

	@Override
	final public void checkBounds(double t, double[] x) {
		if (isOutOfBounds(t, x))
			listener.reactionRateBoundEventOccured(reaction, t, x);
	}

	@Override
	public void init(double t0, double[] x0, double t) {
	}

	@Override
	final public double g(double t, double[] x) {
		double rate = model.computePropensity(reaction, t, x);
		switch (boundType) {
		case LOWER:
			return rate - lowerBound;
		case UPPER:
			return rate - upperBound;
		case BOTH:
			return (rate - lowerBound) * (rate - upperBound);
		case NONE:
		default:
			return 1.0;
		}
	}

	@Override
	public Action eventOccurred(double t, double[] x, boolean increasing) {
		listener.reactionRateBoundEventOccured(reaction, t, x);
		return Action.STOP;
	}

	@Override
	public void resetState(double t, double[] x) {
	}

}
