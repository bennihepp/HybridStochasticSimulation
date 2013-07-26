package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import org.apache.commons.math3.ode.events.EventHandler;
import org.apache.commons.math3.util.FastMath;

import ch.ethz.khammash.hybridstochasticsimulation.math.MathUtilities;
import ch.ethz.khammash.hybridstochasticsimulation.models.StateBoundObserver.BoundType;
import ch.ethz.khammash.hybridstochasticsimulation.networks.AdaptiveMSHRN;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork.ReactionType;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork.SpeciesType;
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
	private double[] optionalState;
	private int parentOptionalStateSize;

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
		numberOfAdapations = 0;
	}

//	public AdaptiveMSHRNModel(AdaptiveMSHRNModel model) {
//		this(model.hrn);
//	}

	@Override
	public void initialize(double t0, double[] x0) {
		hrn.init();
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
			// TODO: Multiply by x oder not?
			double lowerBound1 = x * Math.pow(hrn.getN(), -hrn.getEta());
			double lowerBound2 = hrn.getInverseSpeciesScaleFactor(s) * FastMath.pow(hrn.getN(), hrn.getXi() - hrn.getEta());
			double lowerBound = Math.max(lowerBound1, lowerBound2);
			lowerBound = Math.min(lowerBound, x);
			double upperBound = x * Math.pow(hrn.getN(),  hrn.getEta());
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
		for (int r=0; r < getNumberOfReactions(); r++)
			tmpPropensities[r] = computePropensity(r, t, x);
		computeDerivatives(t, x, tmpxDot);
		hrn.adapt(t, x, tmpxDot, tmpPropensities);
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

	@Override
	public double[] computeOptionalState(double t, double[] x) {
		double[] tmp = super.computeOptionalState(t, x);
		parentOptionalStateSize = tmp.length;
		if (optionalState == null)
			optionalState = new double[getNumberOfSpecies() + getNumberOfReactions() + tmp.length];
		int i = 0;
		for (int j=0; j < tmp.length; j++)
			optionalState[i++] = tmp[j];
		for (int s=0; s < getNumberOfSpecies(); s++)
			optionalState[i++] = computeSpeciesType(s);
		for (int r=0; r < getNumberOfReactions(); r++)
			optionalState[i++] = computeReactionType(r);
		return optionalState;
	}

	public List<Integer> getSpeciesTypeStateIndices() {
		int from = parentOptionalStateSize;
		int to = from + getNumberOfSpecies();
		return MathUtilities.intRangeList(from, to);
	}

	public List<Integer> getReactionTypeStateIndices() {
		int from = parentOptionalStateSize + getNumberOfSpecies();
		int to = from + getNumberOfReactions();
		return MathUtilities.intRangeList(from, to);
	}

	private double computeReactionType(int reaction) {
		ReactionType reactionType = hrn.getReactionType(reaction);
		switch (reactionType) {
		case NONE:
			return 0.0;
		case STOCHASTIC:
			return -(reaction + 1);
		case DETERMINISTIC:
			return +(reaction + 1);
		case EXPLODING:
			return (reaction + 1) + getNumberOfReactions();
		default:
			return Double.NaN;
		}
	}

	private double computeSpeciesType(int species) {
		SpeciesType speciesType = hrn.getSpeciesType(species);
		switch (speciesType) {
		case CONTINUOUS:
			return +(species + 1);
		case DISCRETE:
			return -(species + 1);
		case UNDEFINED:
			return 0;
		default:
			return Double.NaN;
		}
	}

	public void setOptionalEventOccured() {
		hasOptionalEventOccured = true;
	}

}
