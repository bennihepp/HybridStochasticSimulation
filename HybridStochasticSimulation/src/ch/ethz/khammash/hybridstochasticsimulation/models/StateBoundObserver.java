package ch.ethz.khammash.hybridstochasticsimulation.models;

import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPEventObserver;


public class StateBoundObserver implements PDMPEventObserver {

	public static enum BoundType {
		LOWER, UPPER, BOTH,
	}

	private StateBoundEventListener listener;
	private int species;
	private double lowerBound;
	private double upperBound;
	private BoundType boundType;

	public StateBoundObserver(StateBoundEventListener listener, int species, double bound, BoundType boundType) {
		this.listener = listener;
		this.species = species;
		lowerBound = bound;
		upperBound = bound;
		this.boundType = boundType;
	}

	public StateBoundObserver(StateBoundEventListener listener, int species, double min, double max) {
		this.listener = listener;
		this.species = species;
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

	final public boolean isOutOfBounds(double[] x) {
		switch (boundType) {
		case LOWER:
			return x[species] <= lowerBound;
		case UPPER:
			return x[species] >= upperBound;
		case BOTH:
		default:
			return (x[species] <= lowerBound) || (x[species] >= upperBound);
		}
	}

	final public void checkBounds(double t, double[] x) {
		if (isOutOfBounds(x))
			listener.stateBoundEventOccured(species, t, x);
	}

	@Override
	public void init(double t0, double[] x0, double t) {
	}

	@Override
	final public double g(double t, double[] x) {
		switch (boundType) {
		case LOWER:
			return x[species] - lowerBound;
		case UPPER:
			return x[species] - upperBound;
		case BOTH:
		default:
			return (x[species] - lowerBound) * (x[species] - upperBound);
		}
	}

	@Override
	public Action eventOccurred(double t, double[] x, boolean increasing) {
		listener.stateBoundEventOccured(species, t, x);
		return Action.STOP;
	}

	@Override
	public void resetState(double t, double[] x) {
	}

}
