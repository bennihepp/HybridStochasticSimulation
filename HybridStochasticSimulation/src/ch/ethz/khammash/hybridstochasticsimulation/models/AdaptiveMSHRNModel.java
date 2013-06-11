package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.apache.commons.math3.ode.events.EventHandler;

import ch.ethz.khammash.hybridstochasticsimulation.models.StateBoundObserver.BoundType;
import ch.ethz.khammash.hybridstochasticsimulation.networks.AdaptiveMSHRN;


public class AdaptiveMSHRNModel extends PDMPModelAdapter implements StateBoundEventListener {

	private AdaptiveMSHRN hrn;
	private MSHybridReactionNetworkModel hrnModel;
	private List<EventHandler> optionalEventHandlers;
	private boolean optionalEventFlag;
//	private int optionalEventSpeciesIndex;

	public AdaptiveMSHRNModel(AdaptiveMSHRN hrn) {
		super(new MSHybridReactionNetworkModel(hrn));
		this.hrn = hrn;
		this.hrnModel = (MSHybridReactionNetworkModel)getHybridModel();
		optionalEventHandlers = new ArrayList<EventHandler>(this.hrn.getNumberOfSpecies());
		for (int s=0; s < this.hrn.getNumberOfSpecies(); s++)
			optionalEventHandlers.add(new StateBoundObserver(this, s, Double.MAX_VALUE, BoundType.UPPER));
		optionalEventFlag = false;
	}

	@Override
	public void initialize(double t0, double[] x0) {
		adapt(t0, x0);
	}

	private void updateOptionalEventHandlers(double[] x) {
		for (int s=0; s < optionalEventHandlers.size(); s++)
			updateOptionalEventHandler(s, x[s]);
	}

	private void updateOptionalEventHandler(int s, double x) {
		EventHandler eh = optionalEventHandlers.get(s);
		StateBoundObserver seh = (StateBoundObserver)eh;
		if (hrn.getAlpha(s) == 0.0) {
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

	@Override
	public boolean getOptionalEventFlag() {
		return optionalEventFlag;
	}

	@Override
	public void handleOptionalEvent(double t, double[] x) {
		if (optionalEventFlag) {
			adapt(t, x);
			optionalEventFlag = false;
		}
	}

	private void adapt(double t, double[] x) {
		hrn.adapt(x);
		hrnModel.update();
		updateOptionalEventHandlers(x);
	}

	@Override
	public void manualCheckOptionalEvent(double t, double[] x) {
		for (EventHandler eh : optionalEventHandlers) {
			StateBoundObserver seh = (StateBoundObserver)eh;
			seh.checkBounds(t, x);
		}
		handleOptionalEvent(t, x);
	}

	@Override
	public Collection<EventHandler> getOptionalEventHandlers() {
		return optionalEventHandlers;
	}

	@Override
	public boolean isTimeIndependent() {
		return true;
	}

	@Override
	public void stateBoundEventOccured(int species, double t, double[] x) {
		optionalEventFlag = true;
//		optionalEventSpeciesIndex = s;
	}

}
