package ch.ethz.bhepp.hybridstochasticsimulation.models;
//package ch.ethz.khammash.hybridstochasticsimulation.models;
//
//import java.util.ArrayList;
//import java.util.Collections;
//import java.util.HashSet;
//import java.util.LinkedList;
//import java.util.List;
//import java.util.Set;
//
//import org.apache.commons.math3.ode.events.EventHandler;
//import org.apache.commons.math3.util.FastMath;
//
//import ch.ethz.khammash.hybridstochasticsimulation.math.MathUtilities;
//import ch.ethz.khammash.hybridstochasticsimulation.models.StateBoundObserver.BoundType;
//import ch.ethz.khammash.hybridstochasticsimulation.networks.AdaptiveMSHRN;
//import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork.ReactionType;
//import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork.SpeciesType;
//import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPEventObserver;
//import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPEventObserverCollector;
//
//
//public class CustomAdaptiveMSHRNModel extends PDMPMSHRNModel implements StateBoundEventListener {
//
//	public Set<Integer> coupledDiscreteSpecies;
//	public Set<Integer> coupledStochasticReactions;
//	public Set<Integer> uncoupledStochasticReactions;
////	double[] tmpPropVec;
//
//	private AdaptiveMSHRN hrn;
//	private List<PDMPEventObserver> stateBoundObservers;
//	private List<PDMPEventObserver> optionalEventObserverList;
//	private boolean hasOptionalEventOccured;
////	private int optionalEventSpeciesIndex;
//	private double[] tmpPropensities;
//	private double[] tmpxDot;
//	private int numberOfAdapations;
//	private double[] optionalState;
//	private int parentOptionalStateSize;
//
//	public CustomAdaptiveMSHRNModel(AdaptiveMSHRN hrn) {
//		super(hrn);
//		this.hrn = hrn;
//		stateBoundObservers = new ArrayList<PDMPEventObserver>(this.hrn.getNumberOfSpecies());
//		for (int s=0; s < this.hrn.getNumberOfSpecies(); s++)
//			stateBoundObservers.add(new StateBoundObserver(this, s, Double.MAX_VALUE, BoundType.UPPER));
//		hasOptionalEventOccured = false;
//		PDMPEventObserverCollector observerCollector = new PDMPEventObserverCollector(this);
//		for (PDMPEventObserver observer : stateBoundObservers)
//			observerCollector.add(observer);
//		optionalEventObserverList = new LinkedList<PDMPEventObserver>();
//		optionalEventObserverList.add(observerCollector);
//		tmpPropensities = new double[getNumberOfReactions()];
//		tmpxDot = new double[getVectorField().getDimension()];
//		numberOfAdapations = 0;
//
//		coupledDiscreteSpecies = new HashSet<>(hrn.getNumberOfSpecies());
//		coupledStochasticReactions = new HashSet<>(hrn.getNumberOfReactions());
//		uncoupledStochasticReactions = new HashSet<>(hrn.getNumberOfReactions());
////		tmpPropVec = new double[hrn.getNumberOfReactions()];
//	}
//
////	public AdaptiveMSHRNModel(AdaptiveMSHRNModel model) {
////		this(model.hrn);
////	}
//
//	@Override
//	public void initialize(double t0, double[] x0) {
//		hrn.init();
//		adapt(t0, x0);
//	}
//
//	private void updateOptionalEventHandlers(double[] x) {
//		for (int s=0; s < stateBoundObservers.size(); s++)
//			updateOptionalEventHandler(s, x[s]);
//	}
//
//	private void updateOptionalEventHandler(int s, double x) {
//		EventHandler eh = stateBoundObservers.get(s);
//		StateBoundObserver seh = (StateBoundObserver)eh;
//		if (hrn.getAlpha(s) == 0.0) {
//			double upperBound = Math.pow(hrn.getN(), hrn.getXi());
//			upperBound = Math.max(upperBound, x + 1);
//			seh.setUpperBound(upperBound);
//			seh.setBoundType(BoundType.UPPER);
//		} else {
//			// TODO: Multiply by x oder not?
//			// TOOD: Which lower bounds
//			double lowerBound1 = x * Math.pow(hrn.getN(), -hrn.getEta());
//			double lowerBound2 = hrn.getInverseSpeciesScaleFactor(s) * FastMath.pow(hrn.getN(), hrn.getXi() - hrn.getEta());
//			double lowerBound = Math.max(lowerBound1, lowerBound2);
//			lowerBound = Math.min(lowerBound, x);
//			double upperBound = x * Math.pow(hrn.getN(),  hrn.getEta());
//			upperBound = Math.max(upperBound, x);
//			seh.setLowerBound(lowerBound);
//			seh.setUpperBound(upperBound);
//			seh.setBoundType(BoundType.BOTH);
//		}
//	}
//
//	public int getNumberOfAdapations() {
//		return numberOfAdapations;
//	}
//
//	public void resetNumberOfAdapations() {
//		numberOfAdapations = 0;
//	}
//
//	@Override
//	public boolean hasOptionalEventOccured() {
//		return hasOptionalEventOccured;
//	}
//
//	@Override
//	public void handleOptionalEvent(double t, double[] x) {
//		if (hasOptionalEventOccured) {
//			adapt(t, x);
//			hasOptionalEventOccured = false;
//		}
//	}
//
//	private List<Integer> getSpeciesInfluencingReaction(int r) {
//		List<Integer> influcencingSpecies = new ArrayList<>(hrn.getNumberOfSpecies());
//		int[] consumptionStoichiometries = hrn.getConsumptionStoichiometries(r);
//		for (int s=0; s < consumptionStoichiometries.length; s++)
//			if (consumptionStoichiometries[s] > 0)
//				influcencingSpecies.add(s);
//		return influcencingSpecies;
//	}
//
//	private List<Integer> getReactionsChangingSpecies(int species) {
//		List<Integer> changingReactions = new ArrayList<>(hrn.getNumberOfReactions());
//		for (int r=0; r < hrn.getNumberOfReactions(); r++) {
//			if (hrn.getStoichiometry(species, r) != 0)
//				changingReactions.add(r);
//		}
//		return changingReactions;
//	}
//
//	private void adapt(double t, double[] x) {
//		for (int r=0; r < getNumberOfReactions(); r++)
//			tmpPropensities[r] = computePropensity(r, t, x);
//		computeDerivatives(t, x, tmpxDot);
//		hrn.adapt(t, x, tmpxDot, tmpPropensities);
//		update();
//		updateOptionalEventHandlers(x);
//		numberOfAdapations++;
//
//		coupledDiscreteSpecies.clear();
//    	for (int r=0; r < hrn.getNumberOfReactions(); r++) {
//			if (hrn.getReactionType(r) == ReactionType.DETERMINISTIC) {
//				List<Integer> coupledSpeciesList = getSpeciesInfluencingReaction(r);
//				for (int s : coupledSpeciesList) {
//					if (hrn.getSpeciesType(s) == SpeciesType.DISCRETE)
//						coupledDiscreteSpecies.add(s);
//				}
//			}
//    	}
//    	coupledStochasticReactions.clear();
//    	uncoupledStochasticReactions.clear();
//    	for (int s : coupledDiscreteSpecies) {
//    		List<Integer> coupledReactionList = getReactionsChangingSpecies(s);
//    		coupledStochasticReactions.addAll(coupledReactionList);
//    	}
//    	for (int r=0; r < hrn.getNumberOfReactions(); r++)
//    		if (!coupledStochasticReactions.contains(r))
//    			uncoupledStochasticReactions.add(r);
//	}
//
//	@Override
//	public int getDimension() {
//		return getNumberOfSpecies() + 2;
//	}
//
//	@Override
//	public void computeDerivatives(double t, double[] x, double[] xDot) {
////		super.computeDerivativesAndPropensities(t, x, xDot, propVector);
////		xDot[xDot.length - 2] = 0.0;
////		xDot[xDot.length - 1] = 0.0;
////		for (int i = 0; i < propVector.length; i++)
////			xDot[xDot.length - 2] += propVector[i];
////		double propSum = computeDerivativesAndPropensitiesSum(t, x, xDot);
//		super.computeDerivativesAndPropensitiesSum(t, x, xDot);
//		double coupledPropensitiesSum = computeCoupledPropensitiesSum(t, x);
//		double uncoupledPropensitiesSum = computeUncoupledPropensitiesSum(t, x);
//		double propSum = coupledPropensitiesSum + uncoupledPropensitiesSum;
//		xDot[xDot.length - 2] = propSum;
//		xDot[xDot.length - 1] = coupledPropensitiesSum;
//	}
//
//	public double computeCoupledPropensitiesSum(double t, double[] x) {
//    	double coupledPropensitiesSum = 0.0;
//    	for (int r : coupledStochasticReactions)
//    		coupledPropensitiesSum += computePropensity(r, t, x);
//		return coupledPropensitiesSum;
//	}
//
//	public double computeUncoupledPropensitiesSum(double t, double[] x) {
//    	double uncoupledPropensitiesSum = 0.0;
//    	for (int r : uncoupledStochasticReactions)
//    		uncoupledPropensitiesSum += computePropensity(r, t, x);
//		return uncoupledPropensitiesSum;
//	}
//
//	@Override
//	public void checkAndHandleOptionalEvent(double t, double[] x) {
//		for (EventHandler eh : stateBoundObservers) {
//			StateBoundObserver seh = (StateBoundObserver)eh;
//			seh.checkBounds(t, x);
//		}
//		handleOptionalEvent(t, x);
//	}
//
//	@Override
//	public List<PDMPEventObserver> getOptionalEventObservers() {
////		return Collections.unmodifiableList(stateBoundObservers);
//		return Collections.unmodifiableList(optionalEventObserverList);
//	}
//
//	@Override
//	public boolean isTimeIndependent() {
//		return true;
//	}
//
//	@Override
//	public void stateBoundEventOccured(int species, double t, double[] x) {
//		hasOptionalEventOccured = true;
////		optionalEventSpeciesIndex = s;
//	}
//
//	@Override
//	public void stateBoundEventOccured(double t, double[] x) {
//		hasOptionalEventOccured = true;
////		optionalEventSpeciesIndex = -1;
//	}
//
//	@Override
//	public double[] computeOptionalState(double t, double[] x) {
//		double[] tmp = super.computeOptionalState(t, x);
//		parentOptionalStateSize = tmp.length;
//		if (optionalState == null)
//			optionalState = new double[getNumberOfSpecies() + getNumberOfReactions() + tmp.length];
//		int i = 0;
//		for (int j=0; j < tmp.length; j++)
//			optionalState[i++] = tmp[j];
//		for (int s=0; s < getNumberOfSpecies(); s++)
//			optionalState[i++] = computeSpeciesType(s);
//		for (int r=0; r < getNumberOfReactions(); r++)
//			optionalState[i++] = computeReactionType(r);
//		return optionalState;
//	}
//
//	public List<Integer> getSpeciesTypeStateIndices() {
//		int from = parentOptionalStateSize;
//		int to = from + getNumberOfSpecies();
//		return MathUtilities.intRangeList(from, to);
//	}
//
//	public List<Integer> getReactionTypeStateIndices() {
//		int from = parentOptionalStateSize + getNumberOfSpecies();
//		int to = from + getNumberOfReactions();
//		return MathUtilities.intRangeList(from, to);
//	}
//
//	private double computeReactionType(int reaction) {
//		ReactionType reactionType = hrn.getReactionType(reaction);
//		switch (reactionType) {
//		case NONE:
//			return 0.0;
//		case STOCHASTIC:
//			return -(reaction + 1);
//		case DETERMINISTIC:
//			return +(reaction + 1);
//		case EXPLODING:
//			return (reaction + 1) + getNumberOfReactions();
//		default:
//			return Double.NaN;
//		}
//	}
//
//	private double computeSpeciesType(int species) {
//		SpeciesType speciesType = hrn.getSpeciesType(species);
//		switch (speciesType) {
//		case CONTINUOUS:
//			return +(species + 1);
//		case DISCRETE:
//			return -(species + 1);
//		case UNDEFINED:
//			return 0;
//		default:
//			return Double.NaN;
//		}
//	}
//
//	public void setOptionalEventOccured() {
//		hasOptionalEventOccured = true;
//	}
//
//}
