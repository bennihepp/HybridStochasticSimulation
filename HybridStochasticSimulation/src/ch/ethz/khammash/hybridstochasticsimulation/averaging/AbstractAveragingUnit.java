package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import ch.ethz.khammash.hybridstochasticsimulation.graphs.ReactionNetworkGraph;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;
import ch.ethz.khammash.hybridstochasticsimulation.models.AdaptiveMSHRNModel;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork.ReactionType;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;

import com.google.common.base.Predicate;
import com.google.common.collect.Sets;

public abstract class AbstractAveragingUnit implements ModularAveragingUnit {

	protected UnaryBinaryReactionNetwork network;
	protected ReactionNetworkGraph graph;
	private Set<SpeciesVertex> importantSpecies;
	private List<SubnetworkDescription> previousSubnetworksToAverage;
	private SubnetworkEnumerator subnetworkEnumerator;

	public AbstractAveragingUnit(UnaryBinaryReactionNetwork network, Set<SpeciesVertex> importantSpecies) {
		this(network, importantSpecies, null);
	}

	public AbstractAveragingUnit(UnaryBinaryReactionNetwork network, Set<SpeciesVertex> importantSpecies, SubnetworkEnumerator subnetworkEnumerator) {
		this.network = network;
		this.importantSpecies = importantSpecies;
		this.graph = network.getGraph();
		this.previousSubnetworksToAverage = null;
		if (subnetworkEnumerator == null)
			subnetworkEnumerator = new AllSubnetworkEnumerator(network);
		setSubnetworkEnumerator(subnetworkEnumerator);
	}

	protected AbstractAveragingUnit() {}

	protected Iterable<SubnetworkDescription> enumerateSubnetworks() {
		return subnetworkEnumerator;
	}

	protected void copyFrom(AbstractAveragingUnit provider) {
		this.network = provider.network;
		this.importantSpecies = provider.importantSpecies;
		this.graph = this.network.getGraph();
		this.previousSubnetworksToAverage = null;
		this.subnetworkEnumerator = provider.subnetworkEnumerator;
	}

	@Override
	public void reset() {
		previousSubnetworksToAverage = null;
	}

	private class TrivialSubnetworkFilter implements Predicate<SubnetworkDescription> {

		@Override
		public boolean apply(SubnetworkDescription subnetwork) {
			if (subnetwork.getSubnetworkSpecies().size() == network.getNumberOfSpecies() || subnetwork.getSubnetworkSpecies().isEmpty())
				return false;
			for (SpeciesVertex v : subnetwork.getSubnetworkSpecies()) {
				// Skip this subnetwork if it contains any important species
				if (importantSpecies.contains(v))
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

	@Override
	public List<SubnetworkDescription> getSubnetworksToAverageAndResampleState(double t, double[] x, Predicate<SubnetworkDescription> filter) {
		List<SubnetworkDescription> averagingCandidates = findAveragingCandidates(t, x, filter);
		List<SubnetworkDescription> subnetworksToAverage = greedySelectSubnetworksToAverage(averagingCandidates);

		if (t > 400)
			subnetworksToAverage = Collections.emptyList();


		if (previousSubnetworksToAverage != null)
			samplePreviouslyAveragedSubnetworks(t, x, subnetworksToAverage, previousSubnetworksToAverage);

		previousSubnetworksToAverage = subnetworksToAverage;
		return subnetworksToAverage;
	}

	static List<SubnetworkDescription> greedySelectSubnetworksToAverage(List<SubnetworkDescription> averagingCandidates) {
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
		HashSet<SpeciesVertex> speciesToAverage = new HashSet<SpeciesVertex>();
		List<SubnetworkDescription> subnetworksToAverage = new LinkedList<SubnetworkDescription>();
		for (SubnetworkDescription candidate : averagingCandidates)
			if (Sets.intersection(speciesToAverage, candidate.getSubnetworkSpecies()).isEmpty()) {
				speciesToAverage.addAll(candidate.getSubnetworkSpecies());
				subnetworksToAverage.add(candidate);
			}
		// Return the set of subnetworks chosen for averaging
		return subnetworksToAverage;
	}

	public void samplePreviouslyAveragedSubnetworks(double t, double[] x,
			List<SubnetworkDescription> subnetworksToAverage, List<SubnetworkDescription> previousSubnetworksToAverage) {
		// Resample states that have been averaged before but are no longer averaged
		Set<SpeciesVertex> allAveragingSpecies = new HashSet<SpeciesVertex>();
		for (SubnetworkDescription subnetwork : subnetworksToAverage)
			allAveragingSpecies.addAll(subnetwork.getSubnetworkSpecies());
		for (SubnetworkDescription subnetwork : previousSubnetworksToAverage)
			if (!allAveragingSpecies.containsAll(subnetwork.getSubnetworkSpecies()))
				sampleSubnetworkState(t, x, subnetwork);
	}

	public void updateAveraging(AdaptiveMSHRNModel model, double t, double[] x,
			SubnetworkDescription subnetwork) {
		for (int r : subnetwork.getSubnetworkReactions()) {
			model.getNetwork().overrideReactionType(r, ReactionType.NONE);
		}
		for (int r : subnetwork.getSurroundingReactions()) {
			for (int s : subnetwork.getSubnetworkSpeciesIndices()) {
				if (model.getNetwork().getStoichiometry(s, r) != 0) {
					model.getNetwork().overrideReactionType(r, ReactionType.DISCRETE);
					break;
				}
			}
		}
		model.update();
		for (int r : subnetwork.getSubnetworkReactions()) {
			model.setRateParameter(r, 0.0);
		}
		updateSubnetworkState(model, t, x, subnetwork);
	}

	@Override
	public void updateSubnetworkState(AdaptiveMSHRNModel model, double t, double[] x, SubnetworkDescription subnetwork) {
		double[] firstMoments = computeFirstMoments(t, x, subnetwork);
		for (SpeciesVertex vertex : subnetwork.getSubnetworkSpecies()) {
			int s = vertex.getSpecies();
			x[s] = firstMoments[s];
		}
		double[][] secondMoments = computeSecondMoments(t, x, subnetwork);
		for (int r : subnetwork.getSurroundingReactions()) {
			// TODO: ensure pseudo-linearity here?
			int choiceIndex1 = model.getReactionChoiceIndex1(r);
			int choiceIndex2 = model.getReactionChoiceIndex2(r);
			if (subnetwork.containsSpecies(choiceIndex1)
					&& subnetwork.containsSpecies(choiceIndex2)) {
				model.setReactionChoiceIndices(r, -1, -1);
				double rate = model.getRateParameter(r);
				rate *= secondMoments[choiceIndex1][choiceIndex2];
				rate *= model.getNetwork().getInverseSpeciesScaleFactor(choiceIndex1);
				rate *= model.getNetwork().getInverseSpeciesScaleFactor(choiceIndex2);
				model.setRateParameter(r, rate);
			}
			else if (subnetwork.containsSpecies(choiceIndex1)) {
				model.setReactionChoiceIndex1(r, -1);
				double rate = model.getRateParameter(r);
				rate *= firstMoments[choiceIndex1];
				rate *= model.getNetwork().getInverseSpeciesScaleFactor(choiceIndex1);
				model.setRateParameter(r, rate);
			} else if (subnetwork.containsSpecies(choiceIndex2)) {
				model.setReactionChoiceIndex2(r, -1);
				double rate = model.getRateParameter(r);
				rate *= firstMoments[choiceIndex2];
				rate *= model.getNetwork().getInverseSpeciesScaleFactor(choiceIndex2);
				model.setRateParameter(r, rate);
			}
		}
	}

	abstract protected double[] computeFirstMoments(double t, double[] x, SubnetworkDescription subnetwork);

	abstract protected double[][] computeSecondMoments(double t, double[] x, SubnetworkDescription subnetwork);

}
