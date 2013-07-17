package ch.ethz.khammash.hybridstochasticsimulation.networks;

import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import ch.ethz.khammash.hybridstochasticsimulation.sandbox.DependencyGraph;
import ch.ethz.khammash.hybridstochasticsimulation.sandbox.ReactionEdge;
import ch.ethz.khammash.hybridstochasticsimulation.sandbox.ReactionNetworkGraph;
import ch.ethz.khammash.hybridstochasticsimulation.sandbox.SpeciesVertex;

import com.google.common.collect.Sets;

public abstract class AbstractAveragingProvider implements AveragingProvider {

	protected double theta;
	protected UnaryBinaryReactionNetwork network;
	protected ReactionNetworkGraph graph;
	protected Set<SpeciesVertex> importantSpecies;
	protected DependencyGraph depGraph;

	@Override
	public void init(double theta, UnaryBinaryReactionNetwork network, ReactionNetworkGraph graph, Set<SpeciesVertex> importantSpecies) {
		this.theta = theta;
		this.network = network;
		this.graph = graph;
		this.importantSpecies = importantSpecies;
		this.depGraph = new DependencyGraph(graph);
	}

	protected boolean checkAveragingConditions(Set<SpeciesVertex> subnetwork, double[] speciesTimescales) {
		double maxSubnetworkTimescale = computeMaxSubnetworkTimescale(subnetwork, speciesTimescales);
		double minOutsideTimescale = computeMinOutsideTimescale(subnetwork, speciesTimescales);
		return checkAveragingConditions(maxSubnetworkTimescale, minOutsideTimescale);
	}

	protected boolean checkAveragingConditions(double maxSubnetworkTimescale, double minOutsideTimescale) {
//		// Checking for NaNs is implicitly taken care for by the following comparison
//		if (Double.isNaN(maxSubnetworkTimescale) || Double.isNaN(minOutsideTimescale))
//			return false;
		double subnetworkTimescaleRatio = minOutsideTimescale / maxSubnetworkTimescale;
		return subnetworkTimescaleRatio >= theta;
	}

	protected double computeMinOutsideTimescale(Set<SpeciesVertex> subSpecies, double[] speciesTimescales) {
		double minOutsideTimescale = Double.POSITIVE_INFINITY;
		boolean subnetworkIsIsolated = true;
		for (SpeciesVertex v : subSpecies) {
			// Look at the timescale of all outgoing edges that have targets outside of the subnetwork
			for (ReactionEdge edge : graph.outgoingEdgesOf(v)) {
				SpeciesVertex v2 = graph.getEdgeTarget(edge);
				if (subSpecies.contains(v2))
					continue;
				subnetworkIsIsolated = false;
				if (speciesTimescales[v2.getSpecies()] < minOutsideTimescale)
					minOutsideTimescale = speciesTimescales[v2.getSpecies()];
			}
			// Look at the timescale of all incoming edges from outside of the subnetwork where the sources
			// can be influenced by the subnetwork
			for (ReactionEdge edge : graph.incomingEdgesOf(v)) {
				SpeciesVertex v2 = graph.getEdgeSource(edge);
				if (subSpecies.contains(v2))
					continue;
				if (speciesTimescales[v2.getSpecies()] < minOutsideTimescale) {
					for (SpeciesVertex v3 : subSpecies)
						if (depGraph.containsEdge(v3, v2)) {
							minOutsideTimescale = speciesTimescales[v2.getSpecies()];
							break;
						}
				}
			}
		}
		if (subnetworkIsIsolated)
			return Double.NaN;
		return minOutsideTimescale;
	}

	protected double computeMaxSubnetworkTimescale(Set<SpeciesVertex> subSpecies, double[] speciesTimescales) {
		double maxSubnetworkTimescale = 0.0;
		for (SpeciesVertex v : subSpecies) {
			if (speciesTimescales[v.getSpecies()] > maxSubnetworkTimescale)
				maxSubnetworkTimescale = speciesTimescales[v.getSpecies()];
		}
		return maxSubnetworkTimescale;
	}

	protected List<Set<SpeciesVertex>> greedySelectSubnetworksToAverage(List<Set<SpeciesVertex>> averagingCandidates) {
		// Now always choose the candidate subnetworks with the maximum number of species
		// as long as they don't share any species with already chosen subnetworks
		// (this is a simple greedy strategy but should be good enough).
		// Sort candidate subnetworks in decreasing order of their size
		Collections.sort(averagingCandidates, new Comparator<Set<SpeciesVertex>>() {
			@Override
			public int compare(Set<SpeciesVertex> o1, Set<SpeciesVertex> o2) {
				return Integer.compare(o2.size(), o1.size());
			}
		});
		// Make the choices, going from larger to smaller candidate subnetworks.
		HashSet<SpeciesVertex> speciesToAverage = new HashSet<SpeciesVertex>();
		List<Set<SpeciesVertex>> subnetworksToAverage = new LinkedList<Set<SpeciesVertex>>();
		for (Set<SpeciesVertex> candidate : averagingCandidates)
			if (Sets.intersection(speciesToAverage, candidate).isEmpty()) {
				speciesToAverage.addAll(candidate);
				subnetworksToAverage.add(candidate);
			}
		// Return the set of subnetworks chosen for averaging
		return subnetworksToAverage;
	}

	protected void resamplePreviouslyAveragedSubnetworks(double[] x,
			List<Set<SpeciesVertex>> subnetworksToAverage, List<Set<SpeciesVertex>> previousSubnetworksToAverage) {
		// Resample states that have been averaged before but are no longer averaged
		Set<SpeciesVertex> allAveragingSpecies = new HashSet<SpeciesVertex>();
		for (Set<SpeciesVertex> subnetwork : subnetworksToAverage)
			allAveragingSpecies.addAll(subnetwork);
		for (Set<SpeciesVertex> subnetwork : previousSubnetworksToAverage)
			if (!allAveragingSpecies.containsAll(subnetwork))
				resampleFromSteadyStateDistribution(x, subnetwork);
	}

	protected void resampleFromSteadyStateDistribution(double[] x, Set<SpeciesVertex> subnetwork) {}

}
