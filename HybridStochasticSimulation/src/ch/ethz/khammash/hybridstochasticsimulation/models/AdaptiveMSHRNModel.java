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

		public boolean isOutOfBounds(double[] x) {
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
		public double g(double t, double[] x) {
//			double v1 = (x[species] - lowerBound);
//			double v2 = (x[species] - upperBound);
//			double v = - v1 * v2;
//			return v;
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

	private MSHybridReactionNetworkModel hrnModel;
	private List<EventHandler> optionalEventHandlers;
	private boolean optionalEventFlag;
//	private int optionalEventSpeciesIndex;

	public AdaptiveMSHRNModel(MSHybridReactionNetworkModel hrnModel) {
		super(hrnModel);
		this.hrnModel = hrnModel;
		optionalEventHandlers = new ArrayList<EventHandler>(hrnModel.getNumberOfSpecies());
		for (int s=0; s < hrnModel.getNumberOfSpecies(); s++)
//			optionalEventHandlers.add(new StateBoundEventHandler(s, 0, Double.MAX_VALUE));
			optionalEventHandlers.add(new StateBoundEventHandler(s, Double.MAX_VALUE, BoundType.UPPER));
		optionalEventFlag = false;
	}

	@Override
	public void initialize(double t0, double[] x0) {
		for (int s=0; s < optionalEventHandlers.size(); s++) {
			fireOptionalEvent(s);
			handleOptionalEvent(t0, x0);
		}
	}

	private void updateOptionalEventHandler(int s, double x) {
		EventHandler eh = optionalEventHandlers.get(s);
		StateBoundEventHandler seh = (StateBoundEventHandler)eh;
		if (hrnModel.alpha[s] == 0.0) {
			double upperBound = Math.pow(hrnModel.N,  1 - hrnModel.epsilon);
			upperBound = Math.max(upperBound, x);
			seh.setUpperBound(upperBound);
			seh.setBoundType(BoundType.UPPER);
		} else {
			double lowerBound = Math.pow(hrnModel.N,  -hrnModel.epsilon);
			lowerBound = Math.min(lowerBound, x);
			double upperBound = Math.pow(hrnModel.N,  hrnModel.epsilon);
			upperBound = Math.max(upperBound, x);
//			double lowerBound = 0.5;
//			double upperBound = 2.0;
			seh.setLowerBound(lowerBound);
			seh.setUpperBound(upperBound);
			seh.setBoundType(BoundType.BOTH);
		}
	}

	private void fireOptionalEvent(int s) {
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
//			// Update only alpha
//			int s = optionalEventSpeciesIndex;
//			x[s] *= hrnModel.speciesScales[s];
//			double newAlpha;
//			if (x[s] < Math.pow(hrnModel.N,  1 - hrnModel.epsilon))
//				newAlpha = 0.0;
//			else
//				newAlpha = Math.log(x[s]) / Math.log(hrnModel.N);
//			hrnModel.alpha[s] = newAlpha;
//			hrnModel.adaptScales();
//			x[s] *= hrnModel.inverseSpeciesScales[s];
//			updateOptionalEventHandler(s, x[s]);
			// or update all alphas?
			for (int s=0; s < hrnModel.getNumberOfSpecies(); s++) {
				x[s] *= hrnModel.speciesScales[s];
				double newAlpha;
				if (x[s] < Math.pow(hrnModel.N,  1 - hrnModel.epsilon))
					newAlpha = 0.0;
				else
					newAlpha = Math.log(x[s]) / Math.log(hrnModel.N);
				hrnModel.alpha[s] = newAlpha;
			}
			hrnModel.adaptScales();
			for (int s=0; s < hrnModel.getNumberOfSpecies(); s++) {
				x[s] *= hrnModel.inverseSpeciesScales[s];
				if (hrnModel.alpha[s] == 0.0)
					x[s] = Math.round(x[s]);
				updateOptionalEventHandler(s, x[s]);
			}
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
