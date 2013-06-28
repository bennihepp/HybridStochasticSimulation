package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.math3.ode.events.EventHandler;
import org.apache.commons.math3.util.FastMath;

import ch.ethz.khammash.hybridstochasticsimulation.models.StateBoundObserver.BoundType;
import ch.ethz.khammash.hybridstochasticsimulation.networks.AdaptiveMSHRN;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPEventObserver;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPEventObserverCollector;


public class AdaptiveMSHRNModel extends PDMPMSHRNModel implements StateBoundEventListener {

	private AdaptiveMSHRN hrn;
	private List<PDMPEventObserver> stateBoundObservers;
	private List<PDMPEventObserver> optionalEventObserverList;
	private boolean hasOptionalEventOccured;
//	private int optionalEventSpeciesIndex;
	private double[] tmpPropensities;
	private double[] tmpxDot;
	private int numberOfAdapations;

	public AdaptiveMSHRNModel(AdaptiveMSHRN hrn) {
		super(hrn);
		this.hrn = hrn;
		stateBoundObservers = new ArrayList<PDMPEventObserver>(this.hrn.getNumberOfSpecies());
		for (int s=0; s < this.hrn.getNumberOfSpecies(); s++)
			stateBoundObservers.add(new StateBoundObserver(this, s, Double.MAX_VALUE, BoundType.UPPER));
		hasOptionalEventOccured = false;
		PDMPEventObserverCollector observerCollector = new PDMPEventObserverCollector(this);
		for (PDMPEventObserver observer : stateBoundObservers)
			observerCollector.add(observer);
		optionalEventObserverList = new LinkedList<PDMPEventObserver>();
		optionalEventObserverList.add(observerCollector);
		tmpPropensities = new double[getNumberOfReactions()];
		tmpxDot = new double[getNumberOfSpecies()];
		this.numberOfAdapations = 0;
	}

	public AdaptiveMSHRNModel(AdaptiveMSHRNModel model) {
		this(model.hrn);
	}

	@Override
	public void initialize(double t0, double[] x0) {
		adapt(t0, x0);
	}

	private void updateOptionalEventHandlers(double[] x) {
		for (int s=0; s < stateBoundObservers.size(); s++)
			updateOptionalEventHandler(s, x[s]);
	}

	private void updateOptionalEventHandler(int s, double x) {
		EventHandler eh = stateBoundObservers.get(s);
		StateBoundObserver seh = (StateBoundObserver)eh;
		if (hrn.getAlpha(s) == 0.0) {
			double upperBound = Math.pow(hrn.getN(), hrn.getXi());
			upperBound = Math.max(upperBound, x + 1);
			seh.setUpperBound(upperBound);
			seh.setBoundType(BoundType.UPPER);
		} else {
			// TODO: Properly handle switching back to real copy numbers
			// TODO: Multiply by x oder not?
			double lowerBound1 = x * Math.pow(hrn.getN(), -hrn.getEpsilon());
			double lowerBound2 = hrn.getInverseSpeciesScaleFactor(s) * FastMath.pow(hrn.getN(), hrn.getXi() - hrn.getEpsilon());
			double lowerBound = Math.max(lowerBound1, lowerBound2);
			lowerBound = Math.min(lowerBound, x);
			double upperBound = x * Math.pow(hrn.getN(),  hrn.getEpsilon());
			upperBound = Math.max(upperBound, x);
			seh.setLowerBound(lowerBound);
			seh.setUpperBound(upperBound);
			seh.setBoundType(BoundType.BOTH);
		}
	}

	public int getNumberOfAdapations() {
		return numberOfAdapations;
	}

	public void resetNumberOfAdapations() {
		numberOfAdapations = 0;
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
//		for (int r=0; r < getNumberOfReactions(); r++)
//		tmpPropensities[r] = computePropensity(r, t, x);
		computeDerivativesAndPropensities(t, x, tmpxDot, tmpPropensities);
		hrn.adapt(x, tmpxDot, tmpPropensities);
		update();
		updateOptionalEventHandlers(x);
		numberOfAdapations++;
	}

	@Override
	public void checkAndHandleOptionalEvent(double t, double[] x) {
		for (EventHandler eh : stateBoundObservers) {
			StateBoundObserver seh = (StateBoundObserver)eh;
			seh.checkBounds(t, x);
		}
		handleOptionalEvent(t, x);
	}

	@Override
	public List<PDMPEventObserver> getOptionalEventObservers() {
//		return Collections.unmodifiableList(stateBoundObservers);
		return Collections.unmodifiableList(optionalEventObserverList);
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

	@Override
	public void stateBoundEventOccured(double t, double[] x) {
		hasOptionalEventOccured = true;
//		optionalEventSpeciesIndex = -1;
	}

}
