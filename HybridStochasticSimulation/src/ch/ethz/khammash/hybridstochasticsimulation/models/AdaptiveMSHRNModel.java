package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.ode.events.EventHandler;
import org.apache.commons.math3.util.FastMath;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ch.ethz.khammash.hybridstochasticsimulation.averaging.AveragingCandidateFilter;
import ch.ethz.khammash.hybridstochasticsimulation.averaging.AveragingUnit;
import ch.ethz.khammash.hybridstochasticsimulation.averaging.SubnetworkDescription;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;
import ch.ethz.khammash.hybridstochasticsimulation.math.MathUtilities;
import ch.ethz.khammash.hybridstochasticsimulation.models.StateBoundObserver.BoundType;
import ch.ethz.khammash.hybridstochasticsimulation.networks.AdaptiveMSHRN;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork.ReactionType;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork.SpeciesType;
import ch.ethz.khammash.hybridstochasticsimulation.networks.SpeciesTimescaleSeparationSubnetworkFilter;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPEventObserver;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPEventObserverCollector;

import com.google.common.base.Optional;


public class AdaptiveMSHRNModel extends PDMPMSHRNModel implements StateBoundEventListener {

	private static final Logger logger = LoggerFactory.getLogger(AdaptiveMSHRNModel.class);

	public Set<Integer> coupledDiscreteSpecies;
	public Set<Integer> coupledStochasticReactions;
	public Set<Integer> uncoupledStochasticReactions;

	private AdaptiveMSHRN hrn;
	public List<PDMPEventObserver> stateBoundObservers;
	private List<PDMPEventObserver> optionalEventObserverList;
	private boolean hasOptionalEventOccured;
//	private int optionalEventSpeciesIndex;
//	private double[] tmpPropensities;
//	private double[] tmpxDot;
	private int numberOfAdapations;
	private double[] optionalState;
	private int parentOptionalStateSize;

	private Optional<AveragingUnit> averagingUnitOptional;
	private AveragingCandidateFilter subnetworkFilter;

	private List<SubnetworkDescription> subnetworksToAverage;
	private Set<Integer> reactionsToAverage;

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
//		tmpPropensities = new double[getNumberOfReactions()];
//		tmpxDot = new double[getNumberOfSpecies()];
		numberOfAdapations = 0;

		coupledDiscreteSpecies = new HashSet<>(hrn.getNumberOfSpecies());
		coupledStochasticReactions = new HashSet<>(hrn.getNumberOfReactions());
		uncoupledStochasticReactions = new HashSet<>(hrn.getNumberOfReactions());

		averagingUnitOptional = Optional.absent();
		subnetworkFilter = new SpeciesTimescaleSeparationSubnetworkFilter(hrn);
		subnetworksToAverage = Collections.emptyList();
		reactionsToAverage = new HashSet<>();
	}

	// Important: reset() must be called after changing the averaging unit
	public void setAveragingUnit(AveragingUnit averagingUnit) {
		if (averagingUnit == null)
			unsetAveragingUnit();
		else {
			this.averagingUnitOptional = Optional.of(averagingUnit);
			averagingUnit.reset();
		}
	}

	final public void unsetAveragingUnit() {
		this.averagingUnitOptional = Optional.absent();
	}

	public void setSubnetworkFilter(AveragingCandidateFilter subnetworkFilter) {
		this.subnetworkFilter = subnetworkFilter;
	}

	private void averageSubnetworks(double t, double[] x, List<SubnetworkDescription> subnetworksToAverage) {
		// TODO: don't use hrn.getLogMessages()
		if (hrn.getLogMessages()) {
			boolean[] speciesToAverageMask = new boolean[getNumberOfSpecies()];
			for (SubnetworkDescription subnetwork : subnetworksToAverage) {
				for (SpeciesVertex vertex : subnetwork.getSubnetworkSpecies()) {
					speciesToAverageMask[vertex.getSpecies()] = true;
				}
			}
			HashSet<SpeciesVertex> speciesToAverage = new HashSet<SpeciesVertex>();
			for (int s=0; s < getNumberOfSpecies(); s++)
				if (speciesToAverageMask[s])
					speciesToAverage.add(hrn.getGraph().getSpeciesVertex(s));
			if (logger.isInfoEnabled())
				logger.info(" Species to average {}", speciesToAverage);
		}

		reactionsToAverage.clear();
		for (SubnetworkDescription subnetwork : subnetworksToAverage) {
			reactionsToAverage.addAll(subnetwork.getSurroundingReactions());
		}

		for (SubnetworkDescription subnetwork : subnetworksToAverage) {
			updateAveragedReactions(t, x, subnetwork);
		}
	}

	private void updateAveragedReactions(double t, double[] x, SubnetworkDescription subnetwork) {
		AveragingUnit au = averagingUnitOptional.get();
		au.updateAveraging(this, t, x, subnetwork);
//		updateSubnetworkState(t, x, subnetwork);
//		sampleSubnetworkState(t, x, subnetwork);
	}

	private void updateSubnetworkState(double t, double[] x, SubnetworkDescription subnetwork) {
		AveragingUnit au = averagingUnitOptional.get();
		hrn.invalidateReactionTermTypes();
		au.updateAveraging(this, t, x, subnetwork);
		au.sampleSubnetworkState(t, x, subnetwork);
//		au.updateSubnetworkState(this, t, x, subnetwork);
	}

//	private void sampleSubnetworkState(double t, double[] x, SubnetworkDescription subnetwork) {
//		AveragingUnit au = averagingUnitOptional.get();
//		au.sampleSubnetworkState(t, x, subnetwork);
//	}

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
			// TOOD: Which lower bounds
			double lowerBound1 = x * Math.pow(hrn.getN(), -hrn.getEta());
			// TODO:
//			double lowerBound2 = hrn.getInverseSpeciesScaleFactor(s) * FastMath.pow(hrn.getN(), hrn.getXi() - hrn.getEta());
			// Make 0.9 configurable
			double lowerBound2 = hrn.getInverseSpeciesScaleFactor(s) * FastMath.pow(hrn.getN(), 0.9 * hrn.getXi());
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

	public void adapt(double t, double[] x) {
//		for (int r=0; r < getNumberOfReactions(); r++)
//			tmpPropensities[r] = computePropensity(r, t, x);
//		computeDerivatives(t, x, tmpxDot);
//		hrn.adapt(t, x, tmpxDot, tmpPropensities);

		hrn.unscaleStateInPlace(x);

		hrn.adapt(t, x);
		update();

		if (averagingUnitOptional.isPresent()) {
			AveragingUnit averagingUnit = averagingUnitOptional.get();
//			List<Set<SpeciesVertex>> subnetworksToAverage
//				= averagingUnit.getSubnetworksToAverageAndResampleState(t, x, reactionTimescales);
			subnetworksToAverage
				= averagingUnit.getSubnetworksToAverageAndResampleState(t, x, subnetworkFilter.getFilterPredicate(t, x));
			averageSubnetworks(t, x, subnetworksToAverage);

			// TODO: Is this necessary?
			updateCoupling();
		}

		hrn.scaleStateInPlace(x);

		updateOptionalEventHandlers(x);
		numberOfAdapations++;
	}

	public void manualUpdate(double[] propensities) {
		hrn.invalidateReactionTermTypes();
		hrn.updateReactionTypes(propensities);
		update();
		updateCoupling();
	}

	private void updateCoupling() {
    	coupledStochasticReactions.clear();
		coupledDiscreteSpecies.clear();
    	for (int r=0; r < hrn.getNumberOfReactions(); r++) {
			if (hrn.getReactionType(r) == ReactionType.CONTINUOUS) {
				List<Integer> coupledSpeciesList = getSpeciesInfluencingReaction(r);
				for (int s : coupledSpeciesList) {
					if (hrn.getSpeciesType(s) == SpeciesType.DISCRETE)
						coupledDiscreteSpecies.add(s);
				}
			} else if (hrn.getReactionType(r) == ReactionType.DISCRETE) {
				for (int s=0; s < hrn.getNumberOfSpecies(); s++) {
					if (hrn.getSpeciesType(s) == SpeciesType.CONTINUOUS && hrn.getStoichiometry(s, r) != 0)
						coupledStochasticReactions.add(r);
				}
			}
    	}
    	uncoupledStochasticReactions.clear();
    	for (int s : coupledDiscreteSpecies) {
    		List<Integer> coupledReactionList = getReactionsChangingSpecies(s);
    		for (int r : coupledReactionList) {
    			if (hrn.getReactionType(r) == ReactionType.DISCRETE)
    				coupledStochasticReactions.addAll(coupledReactionList);
    		}
    	}
    	for (int r=0; r < hrn.getNumberOfReactions(); r++)
			if (hrn.getReactionType(r) == ReactionType.DISCRETE) {
	    		if (!coupledStochasticReactions.contains(r))
	    			uncoupledStochasticReactions.add(r);
			}
	}

	private List<Integer> getSpeciesInfluencingReaction(int r) {
		List<Integer> influcencingSpecies = new ArrayList<>(hrn.getNumberOfSpecies());
		int[] consumptionStoichiometries = hrn.getConsumptionStoichiometries(r);
		for (int s=0; s < consumptionStoichiometries.length; s++)
			if (consumptionStoichiometries[s] > 0)
				influcencingSpecies.add(s);
		return influcencingSpecies;
	}

	private List<Integer> getReactionsChangingSpecies(int species) {
		List<Integer> changingReactions = new ArrayList<>(hrn.getNumberOfReactions());
		for (int r=0; r < hrn.getNumberOfReactions(); r++) {
			if (hrn.getStoichiometry(species, r) != 0)
				changingReactions.add(r);
		}
		return changingReactions;
	}

	public double computeDeterministicPropensitiesSum(double t, double[] x) {
    	double deterministicPropensitiesSum = 0.0;
    	for (int r : deterministicReactionIndices)
    		deterministicPropensitiesSum += computePropensity(r, t, x);
		return deterministicPropensitiesSum;
	}

	public double computeCoupledPropensitiesSum(double t, double[] x) {
    	double coupledPropensitiesSum = 0.0;
    	for (int r : coupledStochasticReactions)
    		coupledPropensitiesSum += computePropensity(r, t, x);
		return coupledPropensitiesSum;
	}

	public double computeUncoupledPropensitiesSum(double t, double[] x) {
    	double uncoupledPropensitiesSum = 0.0;
    	for (int r : uncoupledStochasticReactions)
    		uncoupledPropensitiesSum += computePropensity(r, t, x);
		return uncoupledPropensitiesSum;
	}

	public void checkOptionalEvent(double t, double[] x) {
		for (EventHandler eh : stateBoundObservers) {
			StateBoundObserver seh = (StateBoundObserver)eh;
			seh.checkBounds(t, x);
		}
	}

	@Override
	public void checkAndHandleOptionalEvent(double t, double[] x) {
		checkOptionalEvent(t, x);
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
		case DISCRETE:
			return -(reaction + 1);
		case CONTINUOUS:
			return +(reaction + 1);
		case UNDEFINED:
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

	public void handleReaction(int reaction, double t, double[] x, double[] deltaX) {
		changeState(reaction, t, deltaX);
		if (isReactionToAverage(reaction)) {
			// TODO: Only resample state for subnetworks that actually changed species
			hrn.unscaleStateInPlace(x);
			for (SubnetworkDescription subnetwork : subnetworksToAverage) {
				updateSubnetworkState(t, x, subnetwork);
			}
			hrn.scaleStateInPlace(x);
		}
	}

	private boolean isReactionToAverage(int reaction) {
		return reactionsToAverage.contains(reaction);
	}

	public double getRateParameter(int reaction) {
		return modelRateParameters[reaction];
	}

	public void setRateParameter(int reaction, double rate) {
		modelRateParameters[reaction] = rate;
	}

	public void setReactionChoiceIndices(int reaction, int choiceIndex1, int choiceIndex2) {
		reactionChoiceIndices1[reaction] = choiceIndex1;
		reactionChoiceIndices2[reaction] = choiceIndex2;
	}

	public int getReactionChoiceIndex1(int reaction) {
		return reactionChoiceIndices1[reaction];
	}

	public int getReactionChoiceIndex2(int reaction) {
		return reactionChoiceIndices2[reaction];
	}

	public void setReactionChoiceIndex1(int reaction, int choiceIndex) {
		reactionChoiceIndices1[reaction] = choiceIndex;
	}

	public void setReactionChoiceIndex2(int reaction, int choiceIndex) {
		reactionChoiceIndices2[reaction] = choiceIndex;
	}

}
