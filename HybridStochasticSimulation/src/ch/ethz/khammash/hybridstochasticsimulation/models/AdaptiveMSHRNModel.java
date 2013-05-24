package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.apache.commons.math3.ode.events.EventHandler;


public class AdaptiveMSHRNModel extends PDMPModelAdapter {

	private class StateBoundEventHandler implements EventHandler {

		private int species;
		private double min;
		private double max;

		public StateBoundEventHandler(int species, double min, double max) {
			this.species = species;
			this.min = min;
			this.max = max;
		}

		@Override
		public Action eventOccurred(double t, double[] x, boolean increasing) {
			fireOptionalEvent(species);
			return Action.STOP;
		}

		@Override
		public double g(double t, double[] x) {
			double v1 = (x[species] - min);
			double v2 = (x[species] - max);
			double v = - v1 * v2;
			return v;
		}

		@Override
		public void init(double t0, double[] x0, double t) {
		}

		@Override
		public void resetState(double t, double[] x) {
		}

		public void setMin(double min) {
			this.min = min;
		}

		public void setMax(double max) {
			this.max = max;
		}

	}

	private MSHybridReactionNetworkModel hrnModel;
	private List<EventHandler> optionalEventHandlers;
	private boolean optionalEventFlag;
	private int optionalEventSpeciesIndex;
	private double optionalEventNewAlpha;

	public AdaptiveMSHRNModel(MSHybridReactionNetworkModel hrnModel) {
		super(hrnModel);
		this.hrnModel = hrnModel;
		optionalEventHandlers = new ArrayList<EventHandler>(hrnModel.getNumberOfSpecies());
		for (int s=0; s < hrnModel.getNumberOfSpecies(); s++) {
			double min = Math.pow(hrnModel.N,  -0.5);
			double max = Math.pow(hrnModel.N,  0.5);
			optionalEventHandlers.add(new StateBoundEventHandler(s, min, max));
		}
//		initOptionalEventHandlers();
		optionalEventFlag = false;
	}

//	private void initOptionalEventHandlers() {
//		for (int s=0; s < optionalEventHandlers.size(); s++)
//			updateOptionalEventHandler(s);
//		optionalEventFlag = false;
//	}

	private void updateAlpha(int s, double alpha) {
		hrnModel.alpha[s] = alpha;
		hrnModel.adaptScales();
	}

	private void updateOptionalEventHandler(int s, double x) {
		EventHandler eh = optionalEventHandlers.get(s);
		StateBoundEventHandler seh = (StateBoundEventHandler)eh;
		double min;
		if (hrnModel.alpha[s] < 0.5)
			min = 0.0;
		else
			min = Math.pow(hrnModel.N,  -0.5);
		min = Math.min(min, x);
		double max = Math.pow(hrnModel.N,  0.5);
		max = Math.max(max, x);
		seh.setMin(min);
		seh.setMax(max);
	}

	private void fireOptionalEvent(int s) {
		optionalEventFlag = true;
		optionalEventSpeciesIndex = s;
	}

	@Override
	public boolean getOptionalEventFlag() {
		return optionalEventFlag;
	}

	@Override
	public void handleOptionalEvent(double t, double[] x) {
		if (optionalEventFlag) {
			int s = optionalEventSpeciesIndex;
			x[s] *= hrnModel.speciesScales[s];
			double newAlpha = Math.log(x[s]) / Math.log(hrnModel.N);
			newAlpha = Math.max(newAlpha, 0.0);
			updateAlpha(s, newAlpha);
			x[s] *= hrnModel.inverseSpeciesScales[s];
			updateOptionalEventHandler(s, x[s]);
			optionalEventFlag = false;
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
