package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.apache.commons.math3.ode.events.EventHandler;



public class AdaptiveMSHRNModel extends PDMPModelAdapter {

	private enum BoundType {
		LOWER, UPPER, BOTH,
	}

	@SuppressWarnings("unused")
	private class StateBoundEventHandler implements EventHandler {

		private int species;
		private double lowerBound;
		private double upperBound;
		private BoundType boundType;

		public StateBoundEventHandler(int species, double bound, BoundType boundType) {
			this.species = species;
			lowerBound = bound;
			upperBound = bound;
			this.boundType = boundType;
		}

		public StateBoundEventHandler(int species, double min, double max) {
			this.species = species;
			this.lowerBound = min;
			this.upperBound = max;
			this.boundType = BoundType.BOTH;
		}

		@Override
		public Action eventOccurred(double t, double[] x, boolean increasing) {
			fireEvent();
			return Action.STOP;
		}

		public void fireEvent() {
			fireOptionalEvent(species);
		}

		final public boolean isOutOfBounds(double[] x) {
			switch (boundType) {
			case LOWER:
				return x[species] < lowerBound;
			case UPPER:
				return x[species] > upperBound;
			case BOTH:
			default:
				return (x[species] < lowerBound) || (x[species] > upperBound);
			}
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
		public void init(double t0, double[] x0, double t) {
		}

		@Override
		public void resetState(double t, double[] x) {
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

	}

	private AdaptiveMSHybridReactionNetwork hrn;
	private MSHybridReactionNetworkModel hrnModel;
	private List<EventHandler> optionalEventHandlers;
	private boolean optionalEventFlag;
//	private int optionalEventSpeciesIndex;

	public AdaptiveMSHRNModel(AdaptiveMSHybridReactionNetwork hrn) {
		super(new MSHybridReactionNetworkModel(hrn));
		this.hrn = hrn;
		this.hrnModel = (MSHybridReactionNetworkModel)getHybridModel();
		optionalEventHandlers = new ArrayList<EventHandler>(this.hrnModel.getStateDimension());
		for (int s=0; s < this.hrnModel.getStateDimension(); s++)
			optionalEventHandlers.add(new StateBoundEventHandler(s, Double.MAX_VALUE, BoundType.UPPER));
		optionalEventFlag = false;
	}

	@Override
	public void initialize(double t0, double[] x0) {
		hrn.reset();
		hrnModel.update();
		updateOptionalEventHandlers(x0);
	}

	private void updateOptionalEventHandlers(double[] x) {
		for (int s=0; s < optionalEventHandlers.size(); s++)
			updateOptionalEventHandler(s, x[s]);
	}

	private void updateOptionalEventHandler(int s, double x) {
		EventHandler eh = optionalEventHandlers.get(s);
		StateBoundEventHandler seh = (StateBoundEventHandler)eh;
		if (hrn.alpha[s] == 0.0) {
			double upperBound = Math.pow(hrn.getN(),  1 - hrn.getEpsilon());
			upperBound = Math.max(upperBound, x);
			seh.setUpperBound(upperBound);
			seh.setBoundType(BoundType.UPPER);
		} else {
			double lowerBound = x * Math.pow(hrn.getN(), -hrn.getEpsilon());
			lowerBound = Math.min(lowerBound, x);
			double upperBound = x * Math.pow(hrn.getN(),  hrn.getEpsilon());
			upperBound = Math.max(upperBound, x);
			seh.setLowerBound(lowerBound);
			seh.setUpperBound(upperBound);
			seh.setBoundType(BoundType.BOTH);
		}
	}

	final public void fireOptionalEvent(int s) {
		optionalEventFlag = true;
//		optionalEventSpeciesIndex = s;
	}

	@Override
	public boolean getOptionalEventFlag() {
		return optionalEventFlag;
	}

	@Override
	public void handleOptionalEvent(double t, double[] x) {
		if (optionalEventFlag) {
			hrn.adapt(x);
			hrnModel.update();
			updateOptionalEventHandlers(x);
			optionalEventFlag = false;
		}
	}

	@Override
	public void manualCheckOptionalEvent(double t, double[] x) {
		for (EventHandler eh : optionalEventHandlers) {
			StateBoundEventHandler seh = (StateBoundEventHandler)eh;
			if (seh.isOutOfBounds(x)) {
				seh.fireEvent();
				handleOptionalEvent(t, x);
				break;
			}
		}
	}

	@Override
	public Collection<EventHandler> getOptionalEventHandlers() {
		return optionalEventHandlers;
	}

	@Override
	public boolean isTimeIndependent() {
		return true;
	}

}
