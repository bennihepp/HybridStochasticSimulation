package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.apache.commons.math3.ode.events.EventHandler;

import ch.ethz.khammash.hybridstochasticsimulation.models.StateBoundObserver.BoundType;
import ch.ethz.khammash.hybridstochasticsimulation.networks.AdaptiveMSHRN;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPEventObserver;


public class AdaptiveMSHRNModel extends PDMPModelAdapter<MSHybridReactionNetworkModel> implements StateBoundEventListener {

	private AdaptiveMSHRN hrn;
	private MSHybridReactionNetworkModel hrnModel;
	private List<PDMPEventObserver> optionalEventObservers;
	private boolean hasOptionalEventOccured;
//	private int optionalEventSpeciesIndex;

	public AdaptiveMSHRNModel(AdaptiveMSHRN hrn) {
		super(new MSHybridReactionNetworkModel(hrn));
		this.hrn = hrn;
		this.hrnModel = (MSHybridReactionNetworkModel)getHybridModel();
		optionalEventObservers = new ArrayList<>(this.hrn.getNumberOfSpecies());
		for (int s=0; s < this.hrn.getNumberOfSpecies(); s++)
			optionalEventObservers.add(new StateBoundObserver(this, s, Double.MAX_VALUE, BoundType.UPPER));
		hasOptionalEventOccured = false;
	}

	@Override
	public void initialize(double t0, double[] x0) {
		adapt(t0, x0);
	}

	private void updateOptionalEventHandlers(double[] x) {
		for (int s=0; s < optionalEventObservers.size(); s++)
			updateOptionalEventHandler(s, x[s]);
	}

	private void updateOptionalEventHandler(int s, double x) {
		EventHandler eh = optionalEventObservers.get(s);
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
	public boolean hasOptionalEventOccured() {
		return hasOptionalEventOccured;
	}

	@Override
	public void handleOptionalEvent(double t, double[] x) {
		if (hasOptionalEventOccured) {
			adapt(t, x);
			hasOptionalEventOccured = false;
		}
	}

	private void adapt(double t, double[] x) {
		hrn.adapt(x);
		hrnModel.update();
		updateOptionalEventHandlers(x);
	}

	@Override
	public void checkForOptionalEvent(double t, double[] x) {
		for (EventHandler eh : optionalEventObservers) {
			StateBoundObserver seh = (StateBoundObserver)eh;
			seh.checkBounds(t, x);
		}
		handleOptionalEvent(t, x);
	}

	@Override
	public Collection<PDMPEventObserver> getOptionalEventObservers() {
		return optionalEventObservers;
	}

	@Override
	public boolean isTimeIndependent() {
		return true;
	}

	@Override
	public void stateBoundEventOccured(int species, double t, double[] x) {
		hasOptionalEventOccured = true;
//		optionalEventSpeciesIndex = s;
	}

}
