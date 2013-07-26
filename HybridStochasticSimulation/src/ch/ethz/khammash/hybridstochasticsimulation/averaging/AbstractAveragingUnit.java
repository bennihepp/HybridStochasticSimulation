package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.util.FastMath;

import ch.ethz.khammash.hybridstochasticsimulation.graphs.ReactionEdge;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.ReactionNetworkGraph;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;
import ch.ethz.khammash.hybridstochasticsimulation.networks.UnaryBinaryReactionNetwork;

import com.google.common.collect.Sets;

public abstract class AbstractAveragingUnit implements AveragingUnit {

	protected double theta;
	protected UnaryBinaryReactionNetwork network;
	protected ReactionNetworkGraph graph;
	protected Set<SpeciesVertex> importantSpecies;
	private List<Set<SpeciesVertex>> previousSubnetworksToAverage;
	private SubnetworksEnumerator subnetworksEnumerator;

	public AbstractAveragingUnit(double theta, UnaryBinaryReactionNetwork network, Set<SpeciesVertex> importantSpecies) {
		init(theta, network, importantSpecies); 
	}

	protected AbstractAveragingUnit() {}

	final private void init(double theta, UnaryBinaryReactionNetwork network, Set<SpeciesVertex> importantSpecies) {
		this.theta = theta;
		this.network = network;
		this.importantSpecies = importantSpecies;
		this.graph = network.getGraph();
		this.previousSubnetworksToAverage = null;
		this.subnetworksEnumerator = new AllSubnetworksEnumerator(graph);
	}

	protected SubnetworksEnumerator getSubnetworksEnumerator() {
		return subnetworksEnumerator;
	}

	protected void copyFrom(AbstractAveragingUnit provider) {
		init(provider.theta, provider.network, provider.importantSpecies);
	}

	@Override
	public void reset() {
		previousSubnetworksToAverage = null;
	}

	@Override
	public void setSubnetworksEnumerator(SubnetworksEnumerator subnetworksEnumerator) {
		this.subnetworksEnumerator = subnetworksEnumerator;
	}

	@Override
	public List<Set<SpeciesVertex>> getSubnetworksToAverageAndResampleState(double t, double[] x, double[] reactionTimescales) {
		List<Set<SpeciesVertex>> averagingCandidates = findAveragingCandidates(t, x, reactionTimescales);
		List<Set<SpeciesVertex>> subnetworksToAverage = greedySelectSubnetworksToAverage(averagingCandidates);

		if (previousSubnetworksToAverage != null)
			resamplePreviouslyAveragedSubnetworks(t, x, subnetworksToAverage, previousSubnetworksToAverage);

		previousSubnetworksToAverage = subnetworksToAverage;
		return subnetworksToAverage;
	}

	protected boolean checkAveragingConditions(Set<SpeciesVertex> subnetworkSpecies, double[] x, double[] reactionTimescales) {
		double maxSubnetworkTimescale = computeMaxSubnetworkTimescale(subnetworkSpecies, x, reactionTimescales);
		double minOutsideTimescale = computeMinBorderTimescale(subnetworkSpecies, x, reactionTimescales);
		boolean result = checkAveragingConditions(maxSubnetworkTimescale, minOutsideTimescale);
		return result;
	}

	protected boolean checkAveragingConditions(double maxSubnetworkTimescale, double minOutsideTimescale) {
//		// Checking for NaNs is implicitly taken care for by the following comparison
//		if (Double.isNaN(maxSubnetworkTimescale) || Double.isNaN(minOutsideTimescale))
//			return false;
		double subnetworkTimescaleRatio = minOutsideTimescale / maxSubnetworkTimescale;
		boolean result = subnetworkTimescaleRatio >= theta;
		return result;
	}

	protected double computeMinBorderTimescale(Set<SpeciesVertex> subnetworkSpecies, double[] x, double[] reactionTimescales) {
		Set<ReactionEdge> borderEdges = new HashSet<ReactionEdge>();
		for (SpeciesVertex v : subnetworkSpecies) {
			for (ReactionEdge e : graph.incomingEdgesOf(v)) {
				if (!subnetworkSpecies.contains(e.getSource()) && !borderEdges.contains(e))
					borderEdges.add(e);
			}
			for (ReactionEdge e : graph.outgoingEdgesOf(v)) {
				if (!subnetworkSpecies.contains(e.getTarget()) && !borderEdges.contains(e))
					borderEdges.add(e);
			}
		}
		double minBorderTimescale = Double.POSITIVE_INFINITY;;
		for (ReactionEdge e : borderEdges) {
			double timescale = reactionTimescales[e.getReaction()];
			double xOutside;
			if (subnetworkSpecies.contains(e.getSource()))
				xOutside = x[e.getTarget().getSpecies()];
			else
				xOutside = x[e.getSource().getSpecies()];
			timescale /= FastMath.max(xOutside, 1.0);
			if (timescale < minBorderTimescale)
				minBorderTimescale = timescale;
		}
		return minBorderTimescale;
	}

	protected double computeMaxSubnetworkTimescale(Set<SpeciesVertex> subnetworkSpecies, double[] x, double[] reactionTimescales) {
		Set<ReactionEdge> subnetworkEdges = new HashSet<ReactionEdge>();
		for (SpeciesVertex v : subnetworkSpecies) {
			for (ReactionEdge e : graph.incomingEdgesOf(v)) {
				if (subnetworkSpecies.contains(e.getSource()) && !subnetworkEdges.contains(e))
					subnetworkEdges.add(e);
			}
			for (ReactionEdge e : graph.outgoingEdgesOf(v)) {
				if (subnetworkSpecies.contains(e.getTarget()) && !subnetworkEdges.contains(e))
					subnetworkEdges.add(e);
			}
		}
		if (subnetworkEdges.size() == 0)
			return Double.POSITIVE_INFINITY;
		double maxSubnetworkTimescale = 0.0;
		for (ReactionEdge e : subnetworkEdges) {
			double timescale = reactionTimescales[e.getReaction()];
			double xSource = x[e.getSource().getSpecies()];
			double xTarget = x[e.getTarget().getSpecies()];
			double xMin = FastMath.min(xSource, xTarget);
			timescale /= FastMath.max(xMin, 1.0);
			if (timescale > maxSubnetworkTimescale)
				maxSubnetworkTimescale = timescale;
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

	public void resamplePreviouslyAveragedSubnetworks(double t, double[] x,
			List<Set<SpeciesVertex>> subnetworksToAverage, List<Set<SpeciesVertex>> previousSubnetworksToAverage) {
		// Resample states that have been averaged before but are no longer averaged
		Set<SpeciesVertex> allAveragingSpecies = new HashSet<SpeciesVertex>();
		for (Set<SpeciesVertex> subnetwork : subnetworksToAverage)
			allAveragingSpecies.addAll(subnetwork);
		for (Set<SpeciesVertex> subnetworkSpecies : previousSubnetworksToAverage)
			if (!allAveragingSpecies.containsAll(subnetworkSpecies))
				resampleFromSteadyStateDistribution(t, x, subnetworkSpecies);
	}

}
