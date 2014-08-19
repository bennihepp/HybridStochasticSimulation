package ch.ethz.bhepp.hybridstochasticsimulation.models;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.ode.events.EventHandler;
import org.apache.commons.math3.util.FastMath;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ch.ethz.bhepp.hybridstochasticsimulation.averaging.AveragingUnit;
import ch.ethz.bhepp.hybridstochasticsimulation.averaging.SubnetworkDescription;
import ch.ethz.bhepp.hybridstochasticsimulation.math.MathUtilities;
import ch.ethz.bhepp.hybridstochasticsimulation.models.StateBoundObserver.BoundType;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.AdaptiveMSHRN;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.MSHybridReactionNetwork.ReactionType;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.MSHybridReactionNetwork.SpeciesType;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.SpeciesTimescaleSeparationFunction;
import ch.ethz.bhepp.hybridstochasticsimulation.simulators.FilterUtilities;
import ch.ethz.bhepp.hybridstochasticsimulation.simulators.PDMPEventObserver;
import ch.ethz.bhepp.hybridstochasticsimulation.simulators.PDMPEventObserverCollector;

import com.google.common.base.Optional;
import com.google.common.base.Predicate;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Sets;
import com.google.common.collect.Sets.SetView;


public class AdaptiveMSHRNModel extends PDMPMSHRNModel implements StateBoundEventListener, ReactionRateBoundEventListener {

	private static final Logger logger = LoggerFactory.getLogger(AdaptiveMSHRNModel.class);

	public Set<Integer> coupledDiscreteSpecies;
	public Set<Integer> coupledStochasticReactions;
	public Set<Integer> uncoupledStochasticReactions;

	private AdaptiveMSHRN hrn;
	private double observationTime;
	public List<PDMPEventObserver> stateBoundObservers;
	private List<PDMPEventObserver> optionalEventObserverList;
	private boolean hasOptionalEventOccured;
//	private int optionalEventSpeciesIndex;
//	private double[] tmpPropensities;
//	private double[] tmpxDot;
	private int numberOfAdaptations;
	private double[] optionalState;
	private int parentOptionalStateSize;

	private Optional<AveragingUnit> averagingUnitOptional;
//	private SubnetworkEnumerator subnetworkEnumerator;
	private List<SubnetworkDescription> subnetworksToAverage;
	private Set<Integer> reactionsToAverage;
	private Predicate<Set<Integer>> reactionSubsetFilter;

	public AdaptiveMSHRNModel(AdaptiveMSHRN hrn, double observationTime) {
		super(hrn);
		this.hrn = hrn;
		this.observationTime = observationTime;
		stateBoundObservers = new ArrayList<PDMPEventObserver>(this.hrn.getNumberOfSpecies());
		for (int s=0; s < this.hrn.getNumberOfSpecies(); s++)
			stateBoundObservers.add(new StateBoundObserver(this, s));
		for (int r=0; r < this.hrn.getNumberOfReactions(); r++)
			stateBoundObservers.add(new ReactionRateBoundObserver(this, this, r));
		hasOptionalEventOccured = false;
		PDMPEventObserverCollector observerCollector = new PDMPEventObserverCollector(this);
		for (PDMPEventObserver observer : stateBoundObservers)
			observerCollector.add(observer);
		optionalEventObserverList = new LinkedList<PDMPEventObserver>();
		optionalEventObserverList.add(observerCollector);
//		tmpPropensities = new double[getNumberOfReactions()];
//		tmpxDot = new double[getNumberOfSpecies()];
		numberOfAdaptations = 0;

		coupledDiscreteSpecies = new HashSet<>(hrn.getNumberOfSpecies());
		coupledStochasticReactions = new HashSet<>(hrn.getNumberOfReactions());
		uncoupledStochasticReactions = new HashSet<>(hrn.getNumberOfReactions());

		averagingUnitOptional = Optional.absent();
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

	// TODO
//	public void setTimescaleSeparationFilter(AveragingCandidateFilter timescaleSeparationFilter) {
//		this.timescaleSeparationFunction = timescaleSeparationFilter;
//	}

	private void averageSubnetworks(double t, double[] x, List<SubnetworkDescription> subnetworksToAverage) {
		// TODO: don't use hrn.getLogMessages()
		if (hrn.getLogMessages() && logger.isInfoEnabled()) {
			if (subnetworksToAverage.isEmpty())
				logger.info(" No subnetworks to average");
			else {
				logger.info(" Subnetworks to average:");
				for (SubnetworkDescription subnetwork : subnetworksToAverage) {
					logger.info("  {}", subnetwork);
				}
			}
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
//		hrn.invalidateReactionTermTypes();
		au.updateSubnetworkState(this, t, x, subnetwork);
//		au.sampleSubnetworkState(t, x, subnetwork);
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
//		for (int s=0; s < stateBoundObservers.size(); s++)
		for (int s=0; s < getNumberOfSpecies(); s++)
			updateStateBoundObserver(s, x[s]);
		for (int r=0; r < getNumberOfReactions(); r++) {
			// FIXME
			double t = 0.0;
			double rate = computePropensity(r, t, x);
			updateReactionRateBoundObserver(r, rate);
		}
	}

	private void updateReactionRateBoundObserver(int r, double rate) {
		EventHandler eh = stateBoundObservers.get(r + getNumberOfSpecies());
		ReactionRateBoundObserver reh = (ReactionRateBoundObserver)eh;
		if (hrn.getReactionType(r) == ReactionType.CONTINUOUS && rate != 0.0) {
    		double lowerBound = rate * FastMath.pow(hrn.getN(), -hrn.getEta());
    		double upperBound = rate * FastMath.pow(hrn.getN(),  hrn.getEta());
    		reh.setLowerBound(lowerBound);
    		reh.setUpperBound(upperBound);
    		reh.setBoundType(ReactionRateBoundObserver.BoundType.BOTH);
		} else {
		    reh.setBoundType(ReactionRateBoundObserver.BoundType.NONE);
		}
	}

	private void updateStateBoundObserver(int s, double x) {
		EventHandler eh = stateBoundObservers.get(s);
		StateBoundObserver seh = (StateBoundObserver)eh;
		if (hrn.getAlpha(s) == 0.0) {
			double upperBound = FastMath.pow(hrn.getN(), hrn.getMu());
			upperBound = FastMath.max(upperBound, x + 1);
//	        logger.info("species {}: x={}, upperBound={}", s, x, upperBound);
			seh.setUpperBound(upperBound);
			seh.setBoundType(BoundType.UPPER);
		} else {
			double lowerBound1 = x * FastMath.pow(hrn.getN(), -hrn.getEta());
			double lowerBound2 = hrn.getInverseSpeciesScaleFactor(s) * FastMath.pow(hrn.getN(), hrn.getMu() - hrn.getEta());
			double lowerBound = FastMath.max(lowerBound1, lowerBound2);
			lowerBound = FastMath.min(lowerBound, x);
			double upperBound = x * FastMath.pow(hrn.getN(), hrn.getEta());
			upperBound = FastMath.max(upperBound, x);
			seh.setLowerBound(lowerBound);
			seh.setUpperBound(upperBound);
			seh.setBoundType(BoundType.BOTH);
		}
	}

	public int getNumberOfAdapations() {
		return numberOfAdaptations;
	}

	public void resetNumberOfAdapations() {
		numberOfAdaptations = 0;
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

	class SubnetworkEnumerator {

		public Iterable<SubnetworkDescription> getIterable(final Set<Integer> species) {
			return new Iterable<SubnetworkDescription>() {

				@Override
				public Iterator<SubnetworkDescription> iterator() {
					Set<Integer> subnetworkReactions = new HashSet<>();
					for (int s : species) {
						subnetworkReactions.addAll(hrn.getInvolvingReactions(s));
					}
//					Set<Integer> subnetworkReactions = new HashSet<>();
//					for (int reaction : allReactions) {
//						for (int s : species) {
//							if (hrn.getStoichiometry(s, reaction) != 0) {
//								subnetworkReactions.add(reaction);
//								break;
//							}
//						}
//					}
					Set<Set<Integer>> reactionsPowerset = Sets.powerSet(subnetworkReactions);
					final Iterator<Set<Integer>> baseIterator = reactionsPowerset.iterator();
					Iterator<SubnetworkDescription> iterator = new Iterator<SubnetworkDescription>() {

						@Override
						public boolean hasNext() {
							return baseIterator.hasNext();
						}

						@Override
						public SubnetworkDescription next() {
							Set<Integer> reactions = baseIterator.next();
							return new SubnetworkDescription(reactions, hrn);
						}

						@Override
						public void remove() {
							baseIterator.remove();
						}

					};
					return iterator;
				}

			};
		}
	}

	class SpeciesSubsetEnumerator implements Iterable<Set<Integer>> {

		@Override
		public Iterator<Set<Integer>> iterator() {
			Set<Integer> allSpecies = new HashSet<>();
			for (int s=0; s < hrn.getNumberOfSpecies(); s++)
				allSpecies.add(s);
			Set<Set<Integer>> speciesPowerset = Sets.powerSet(allSpecies);
			Iterator<Set<Integer>> iterator = speciesPowerset.iterator();
			return iterator;
		}

	}

	class SubnetworkInformation extends SubnetworkDescription {

		public double subnetworkTau;
		public double surroundingTau;
		public double deltaTau;

		public SubnetworkInformation(SubnetworkDescription subnetwork, double subnetworkTau, double surroundingTau) {
			super(subnetwork);
			this.subnetworkTau = subnetworkTau;
			this.surroundingTau = surroundingTau;
			deltaTau = surroundingTau - subnetworkTau;
		}

		public String toString() {
			return String.format("%s, tau=[%f, %f], deltaTau=%f", super.toString(), subnetworkTau, surroundingTau, deltaTau);
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

		// Averaging happens here
		if (averagingUnitOptional.isPresent()) {
			AveragingUnit au = averagingUnitOptional.get();

			final double[] reactionTimescales = hrn.computeReactionTimescales(t, x);
			List<Integer> allReactions = MathUtilities.intRangeList(0, hrn.getNumberOfReactions());
			Collections.sort(allReactions, new Comparator<Integer>() {

				@Override
				public int compare(Integer o1, Integer o2) {
					return Double.compare(reactionTimescales[o2], reactionTimescales[o1]);
				}

			});

			// TODO: comment
			List<SubnetworkInformation> subnetworkCandidates = new LinkedList<>();
			for (int i=0; i < hrn.getNumberOfReactions() - 1; i++) {
				double deltaTau = reactionTimescales[allReactions.get(i)] - reactionTimescales[allReactions.get(i + 1)];
				if (deltaTau > hrn.getTheta()) {
					Set<Integer> reactionSet = ImmutableSet.copyOf(allReactions.subList(0, i + 1));
					Set<Set<Integer>> reactionPowerSet = Sets.powerSet(reactionSet);
					for (Set<Integer> subnetworkReactions : FilterUtilities.filter(reactionPowerSet, reactionSubsetFilter)) {
						SubnetworkDescription subnetwork = new SubnetworkDescription(subnetworkReactions, hrn);
						if (!Sets.intersection(subnetwork.getSurroundingReactions(), reactionSet).isEmpty())
							continue;
						if (au.getSubnetworkFilter().apply(subnetwork)) {
							double subnetworkTau = SpeciesTimescaleSeparationFunction.computeMinTimescale(subnetwork.getSubnetworkReactions(), reactionTimescales);
							double surroundingTau = SpeciesTimescaleSeparationFunction.computeMaxTimescale(subnetwork.getSurroundingReactions(), reactionTimescales);
							double observationTau = hrn.computeObservationTimescale(observationTime);
							surroundingTau = FastMath.max(surroundingTau, observationTau);
							SubnetworkInformation subnetworkInfo = new SubnetworkInformation(subnetwork, subnetworkTau, surroundingTau);
							subnetworkCandidates.add(subnetworkInfo);
						}
					}
				}
			}

			List<SubnetworkDescription> previousSubnetworksToAverage = subnetworksToAverage;
			subnetworksToAverage = greedySelectSubnetworksToAverage(subnetworkCandidates);
			samplePreviouslyAveragedSubnetworks(t, x, au, subnetworksToAverage, previousSubnetworksToAverage);
			au.reset();
			averageSubnetworks(t, x, subnetworksToAverage);

			// TODO: Is this necessary?
			updateCoupling();
		}

		hrn.scaleStateInPlace(x);

		updateOptionalEventHandlers(x);
		numberOfAdaptations++;
	}

	private void samplePreviouslyAveragedSubnetworks(double t, double[] x,
			AveragingUnit au,
			List<SubnetworkDescription> subnetworksToAverage,
			List<SubnetworkDescription> previousSubnetworksToAverage) {
		// Resample states that have been averaged before but are no longer averaged
		Set<Integer> allAveragingSpecies = new HashSet<>();
		for (SubnetworkDescription subnetwork : subnetworksToAverage)
			allAveragingSpecies.addAll(subnetwork.getSubnetworkSpecies());
		for (SubnetworkDescription subnetwork : previousSubnetworksToAverage)
			if (!allAveragingSpecies.containsAll(subnetwork.getSubnetworkSpecies())) {
				au.sampleSubnetworkState(t, x, subnetwork);
			}
	}

	private List<SubnetworkDescription> greedySelectSubnetworksToAverage(List<SubnetworkInformation> averagingCandidates) {
		// Now always choose the candidate subnetworks with the maximum number of species
		// as long as they don't share any species with already chosen subnetworks
		// (this is a simple greedy strategy but should be good enough).
		// Sort candidate subnetworks in decreasing order of their size
		Collections.sort(averagingCandidates, new Comparator<SubnetworkInformation>() {

			@Override
			public int compare(SubnetworkInformation o1, SubnetworkInformation o2) {
				return Double.compare(o2.subnetworkTau, o1.subnetworkTau);
//				return Integer.compare(o2.getSubnetworkSpecies().size(), o1.getSubnetworkSpecies().size());
			}

		});
		// Make the choices, going from larger to smaller candidate subnetworks.
		HashSet<Integer> speciesToAverage = new HashSet<>();
		HashSet<Integer> reactionsToAverage = new HashSet<>();
		List<SubnetworkDescription> subnetworksToAverage = new LinkedList<SubnetworkDescription>();
		for (SubnetworkInformation candidate : averagingCandidates) {
			SetView<Integer> speciesIntersection = Sets.intersection(speciesToAverage, candidate.getSubnetworkSpecies());
			if (!speciesIntersection.isEmpty())
				continue;
			SetView<Integer> reactionIntersection = Sets.intersection(reactionsToAverage, candidate.getSubnetworkReactions());
			if (!reactionIntersection.isEmpty())
				continue;
			speciesToAverage.addAll(candidate.getSubnetworkSpecies());
			reactionsToAverage.addAll(candidate.getSubnetworkReactions());
			subnetworksToAverage.add(candidate);
		}
		// Return the set of subnetworks chosen for averaging
		return subnetworksToAverage;
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
			OptionalEventObserver obs = (OptionalEventObserver)eh;
			obs.checkBounds(t, x);
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
    public void reactionRateBoundEventOccured(double t, double[] x) {
        hasOptionalEventOccured = true;
//      optionalEventSpeciesIndex = -1;
    }

    @Override
    public void reactionRateBoundEventOccured(int species, double t, double[] x) {
        hasOptionalEventOccured = true;
//      optionalEventSpeciesIndex = s;
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

	public void handleReaction(int reaction, double t, double[] x) {
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

	public void setReactionSubsetFilter(Predicate<Set<Integer>> reactionSubsetFilter) {
		this.reactionSubsetFilter = reactionSubsetFilter;
	}

}
