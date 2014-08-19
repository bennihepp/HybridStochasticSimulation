package ch.ethz.bhepp.hybridstochasticsimulation.averaging;

import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ch.ethz.bhepp.hybridstochasticsimulation.graphs.ReactionNetworkGraph;
import ch.ethz.bhepp.hybridstochasticsimulation.graphs.SpeciesVertex;
import ch.ethz.bhepp.hybridstochasticsimulation.models.AdaptiveMSHRNModel;
import ch.ethz.bhepp.hybridstochasticsimulation.models.UnaryBinaryDeterministicModel;
import ch.ethz.bhepp.hybridstochasticsimulation.models.UnaryBinaryStochasticModel;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.DefaultUnaryBinaryReactionNetwork;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.MSHybridReactionNetwork.ReactionType;
import ch.ethz.bhepp.hybridstochasticsimulation.networks.MassActionReactionNetwork;

import com.google.common.base.Predicate;
import com.google.common.collect.Sets;

public abstract class AbstractAveragingUnit implements ModularAveragingUnit {

	protected MassActionReactionNetwork network;
	protected ReactionNetworkGraph graph;
	private Set<Integer> importantSpecies;
//	private List<SubnetworkDescription> previousSubnetworksToAverage;
	private SubnetworkEnumerator subnetworkEnumerator;
	protected Map<SubnetworkDescription, SubnetworkInformation> subnetworkInformationMap;
	protected int[] reactionChoiceIndices1;
	protected int[] reactionChoiceIndices2;
	protected double[] rateParameters;

	public AbstractAveragingUnit(MassActionReactionNetwork network, Set<Integer> importantSpecies) {
		this(network, importantSpecies, null);
	}

	public AbstractAveragingUnit(MassActionReactionNetwork network, Set<Integer> importantSpecies, SubnetworkEnumerator subnetworkEnumerator) {
		this.network = network;
		this.importantSpecies = importantSpecies;
		this.graph = network.getGraph();
//		this.previousSubnetworksToAverage = null;
		if (subnetworkEnumerator == null)
			subnetworkEnumerator = new AllSubnetworkEnumerator(network);
		setSubnetworkEnumerator(subnetworkEnumerator);
	}

	protected AbstractAveragingUnit() {}

	protected Iterable<SubnetworkDescription> enumerateSubnetworks() {
		return subnetworkEnumerator;
	}

	protected void copyFrom(AbstractAveragingUnit other) {
		this.network = other.network;
		this.importantSpecies = other.importantSpecies;
		this.graph = this.network.getGraph();
//		this.previousSubnetworksToAverage = null;
		this.subnetworkEnumerator = other.subnetworkEnumerator;
		this.subnetworkInformationMap = other.subnetworkInformationMap;
	}

	@Override
	public void reset() {
//		previousSubnetworksToAverage = null;
		this.subnetworkInformationMap = new HashMap<>();
	}

	private class TrivialSubnetworkFilter implements Predicate<SubnetworkDescription> {

		@Override
		public boolean apply(SubnetworkDescription subnetwork) {
			if (subnetwork.getSubnetworkSpecies().size() == network.getNumberOfSpecies() || subnetwork.getSubnetworkSpecies().isEmpty())
				return false;
			for (int s : subnetwork.getSubnetworkSpecies()) {
				// Skip this subnetwork if it contains any important species
				if (importantSpecies.contains(s))
					return false;
			}
			return true;
		}

	}

	@Override
	final public void setSubnetworkEnumerator(SubnetworkEnumerator subnetworkEnumerator) {
		TrivialSubnetworkFilter filter = new TrivialSubnetworkFilter();
		this.subnetworkEnumerator = new FilteredSubnetworkEnumerator(subnetworkEnumerator, filter);
		reset();
	}

//	@Override
//	public List<SubnetworkDescription> getSubnetworksToAverageAndResampleState(double t, double[] x, Predicate<SubnetworkDescription> filter) {
//		List<SubnetworkDescription> averagingCandidates = findAveragingCandidates(t, x, filter);
//		List<SubnetworkDescription> subnetworksToAverage = greedySelectSubnetworksToAverage(averagingCandidates);
//
//		if (t > 400)
//			subnetworksToAverage = Collections.emptyList();
//
//
//		if (previousSubnetworksToAverage != null)
//			samplePreviouslyAveragedSubnetworks(t, x, subnetworksToAverage, previousSubnetworksToAverage);
//
//		previousSubnetworksToAverage = subnetworksToAverage;
//		return subnetworksToAverage;
//	}

	public static List<SubnetworkDescription> greedySelectSubnetworksToAverage(List<SubnetworkDescription> averagingCandidates) {
		// Now always choose the candidate subnetworks with the maximum number of species
		// as long as they don't share any species with already chosen subnetworks
		// (this is a simple greedy strategy but should be good enough).
		// Sort candidate subnetworks in decreasing order of their size
		Collections.sort(averagingCandidates, new Comparator<SubnetworkDescription>() {
			@Override
			public int compare(SubnetworkDescription o1, SubnetworkDescription o2) {
				return Integer.compare(o2.getSubnetworkSpecies().size(), o1.getSubnetworkSpecies().size());
			}
		});
		// Make the choices, going from larger to smaller candidate subnetworks.
		HashSet<Integer> speciesToAverage = new HashSet<>();
		List<SubnetworkDescription> subnetworksToAverage = new LinkedList<SubnetworkDescription>();
		for (SubnetworkDescription candidate : averagingCandidates)
			if (Sets.intersection(speciesToAverage, candidate.getSubnetworkSpecies()).isEmpty()) {
				speciesToAverage.addAll(candidate.getSubnetworkSpecies());
				subnetworksToAverage.add(candidate);
			}
		// Return the set of subnetworks chosen for averaging
		return subnetworksToAverage;
	}

//	public void samplePreviouslyAveragedSubnetworks(double t, double[] x,
//			List<SubnetworkDescription> subnetworksToAverage, List<SubnetworkDescription> previousSubnetworksToAverage) {
//		// Resample states that have been averaged before but are no longer averaged
//		Set<Integer> allAveragingSpecies = new HashSet<>();
//		for (SubnetworkDescription subnetwork : subnetworksToAverage)
//			allAveragingSpecies.addAll(subnetwork.getSubnetworkSpecies());
//		for (SubnetworkDescription subnetwork : previousSubnetworksToAverage)
//			if (!allAveragingSpecies.containsAll(subnetwork.getSubnetworkSpecies())) {
//				SubnetworkInformation subnetworkInfo = subnetworkInformationMap.get(subnetwork);
//				sampleSubnetworkState(t, x, subnetworkInfo);
//			}
//	}

	@Override
	public void sampleSubnetworkState(double t, double[] x, SubnetworkDescription subnetwork) {
		SubnetworkInformation subnetworkInfo = subnetworkInformationMap.get(subnetwork);
		sampleSubnetworkState(t, x, subnetworkInfo);
	}

	abstract protected void sampleSubnetworkState(double t, double[] x, SubnetworkInformation subnetworkInfo);

	protected MassActionReactionNetwork createSubReactionNetwork(MassActionReactionNetwork network, Set<Integer> species, Set<Integer> reactions) {
		DefaultUnaryBinaryReactionNetwork subnetwork = DefaultUnaryBinaryReactionNetwork.createSubnetwork(network, species, reactions);
		return subnetwork;
	}

	protected SubnetworkInformation createSubnetworkInformation(SubnetworkDescription subnetworkDescr) {
		SubnetworkInformation subnetworkInfo = new SubnetworkInformation(subnetworkDescr);
		return subnetworkInfo;
	}

	protected SubnetworkInformation computeSubnetworkInformation(SubnetworkInformation subnetworkInfo) {
		MassActionReactionNetwork subnetwork = createSubReactionNetwork(network, subnetworkInfo.getSubnetworkSpecies(), subnetworkInfo.getSubnetworkReactions());
		subnetworkInfo.setReactionNetwork(subnetwork);

		UnaryBinaryDeterministicModel deterministicModel = new UnaryBinaryDeterministicModel(subnetwork);
		subnetworkInfo.setDeterministicModel(deterministicModel);

		UnaryBinaryStochasticModel stochasticModel = new UnaryBinaryStochasticModel(subnetwork);
		subnetworkInfo.setStochasticModel(stochasticModel);

		List<SpeciesConservationRelation> conservedSpeciesRelations = SpeciesConservationRelation.computeSpeciesConservationRelations(subnetwork);
		subnetworkInfo.setConservedSpeciesRelations(conservedSpeciesRelations);
		Set<Integer> unconservedSpeciesSet = new HashSet<>(subnetworkInfo.getSubnetworkSpecies());
		for (SpeciesConservationRelation relation : subnetworkInfo.getConservedSpeciesRelations())
			unconservedSpeciesSet.removeAll(relation.getConservedSpeciesList());
		subnetworkInfo.setUnconservedSpecies(unconservedSpeciesSet);

		Map<Integer, Integer> subnetworkIndexMap = new HashMap<>();
		for (int i=0; i < subnetwork.getNumberOfSpecies(); i++) {
			SpeciesVertex v = subnetwork.getGraph().getSpeciesVertex(i);
			subnetworkIndexMap.put(v.getSpecies(), i);
		}
		subnetworkInfo.setIndexMap(subnetworkIndexMap);

		return subnetworkInfo;
	}

	public void storeSubnetworkInformation(SubnetworkInformation subnetworkInfo) {
		subnetworkInformationMap.put(subnetworkInfo, subnetworkInfo);
	}

	public void computeAndStoreSubnetworkInformation(SubnetworkDescription subnetworkDescr) {
		SubnetworkInformation subnetworkInfo = createSubnetworkInformation(subnetworkDescr);
		computeSubnetworkInformation(subnetworkInfo);
		storeSubnetworkInformation(subnetworkInfo);
	}

	@Override
	public void updateAveraging(AdaptiveMSHRNModel model, double t, double[] x, SubnetworkDescription subnetworkDescr) {
		computeAndStoreSubnetworkInformation(subnetworkDescr);
		for (int r : subnetworkDescr.getSubnetworkReactions()) {
			model.getNetwork().overrideReactionType(r, ReactionType.NONE);
		}
		for (int r : subnetworkDescr.getSurroundingReactions()) {
			for (int s : subnetworkDescr.getSubnetworkSpecies()) {
				if (model.getNetwork().getStoichiometry(s, r) != 0) {
					model.getNetwork().overrideReactionType(r, ReactionType.DISCRETE);
					break;
				}
			}
		}
		reactionChoiceIndices1 = new int[model.getNumberOfReactions()];
		reactionChoiceIndices2 = new int[model.getNumberOfReactions()];
		rateParameters = new double[model.getNumberOfReactions()];
		for (int r=0; r < model.getNumberOfReactions(); r++) {
			reactionChoiceIndices1[r] = model.getReactionChoiceIndex1(r);
			reactionChoiceIndices2[r] = model.getReactionChoiceIndex2(r);
			rateParameters[r] = model.getRateParameter(r);
		}
		model.update();
		for (int r : subnetworkDescr.getSubnetworkReactions()) {
			model.setRateParameter(r, 0.0);
		}
		updateSubnetworkState(model, t, x, subnetworkDescr);
	}

	@Override
	public void updateSubnetworkState(AdaptiveMSHRNModel model, double t, double[] x, SubnetworkDescription subnetwork) {
		SubnetworkInformation subnetworkInfo = subnetworkInformationMap.get(subnetwork);
		Map<Integer, Integer> indexMap = subnetworkInfo.getIndexMap();

		double[] firstMoments = computeStationaryFirstMoments(t, x, subnetworkInfo);
//		for (int s : subnetwork.getSubnetworkSpecies()) {
//			x[s] = firstMoments[indexMap.get(s)];
//		}
		// TODO: only compute second moments if necessary
		double[][] secondMoments = computeStationarySecondMoments(t, x, subnetworkInfo);

		for (int r : subnetwork.getSurroundingReactions()) {
			// TODO: ensure pseudo-linearity here?
			int choiceIndex1 = reactionChoiceIndices1[r];
			int choiceIndex2 = reactionChoiceIndices2[r];
			double rate = rateParameters[r];
			if (rate < 0 || Double.isNaN(rate))
				throw new RuntimeException();
			if (subnetwork.containsSpecies(choiceIndex1)
					&& subnetwork.containsSpecies(choiceIndex2)) {
				model.setReactionChoiceIndices(r, -1, -1);
				rate *= secondMoments[indexMap.get(choiceIndex1)][indexMap.get(choiceIndex2)];
				rate *= model.getNetwork().getInverseSpeciesScaleFactor(choiceIndex1);
				rate *= model.getNetwork().getInverseSpeciesScaleFactor(choiceIndex2);
				model.setRateParameter(r, rate);
			}
			else if (subnetwork.containsSpecies(choiceIndex1)) {
				if (choiceIndex2 != -1)
					model.setReactionChoiceIndices(r, choiceIndex2, -1);
				else
					model.setReactionChoiceIndex1(r, -1);
				rate *= firstMoments[indexMap.get(choiceIndex1)];
				rate *= model.getNetwork().getInverseSpeciesScaleFactor(choiceIndex1);
				model.setRateParameter(r, rate);
			} else if (subnetwork.containsSpecies(choiceIndex2)) {
				model.setReactionChoiceIndex2(r, -1);
				rate *= firstMoments[indexMap.get(choiceIndex2)];
				rate *= model.getNetwork().getInverseSpeciesScaleFactor(choiceIndex2);
				model.setRateParameter(r, rate);
			}
			if (rate < 0 || Double.isNaN(rate))
				throw new RuntimeException();
		}
		sampleSubnetworkState(t, x, subnetworkInfo);
	}

	abstract protected double[] computeStationaryFirstMoments(double t, double[] x, SubnetworkInformation subnetworkInfo);

	abstract protected double[][] computeStationarySecondMoments(double t, double[] x, SubnetworkInformation subnetworkInfo);

}
