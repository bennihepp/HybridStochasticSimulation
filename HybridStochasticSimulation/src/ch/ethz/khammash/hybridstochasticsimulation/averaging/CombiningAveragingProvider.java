package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ch.ethz.khammash.hybridstochasticsimulation.sandbox.SpeciesVertex;

public class CombiningAveragingProvider extends AbstractAveragingProvider {

	private List<AveragingProvider> averagingProviders;
	private List<Set<SpeciesVertex>> previousSubnetworksToAverage;
	private Map<Set<SpeciesVertex>, AveragingProvider> previousSubnetworkToAveragingProviderMap;

	public CombiningAveragingProvider() {
		averagingProviders = new LinkedList<>();
	}

	public void addAveragingProvider(AveragingProvider ap) {
		averagingProviders.add(ap);
	}

	public void removeAveragingProvider(AveragingProvider ap) {
		averagingProviders.remove(ap);
	}

	public void clearAveragingProviders() {
		averagingProviders.clear();
	}

	@Override
	public List<Set<SpeciesVertex>> getSubnetworksToAverageAndResampleState(double t, double[] x, double[] reactionTimescales) {
		List<Set<SpeciesVertex>> averagingCandidates = findAveragingCandidates(reactionTimescales);
		List<Set<SpeciesVertex>> subnetworksToAverage = greedySelectSubnetworksToAverage(averagingCandidates);

		if (previousSubnetworksToAverage != null)
			resamplePreviouslyAveragedSubnetworks(x, subnetworksToAverage, previousSubnetworksToAverage);

		previousSubnetworksToAverage = subnetworksToAverage;
		return subnetworksToAverage;
	}

	@Override
	public List<Set<SpeciesVertex>> findAveragingCandidates(double[] reactionTimescales) {
		previousSubnetworkToAveragingProviderMap = new HashMap<>();
		List<Set<SpeciesVertex>> allCandidates = new ArrayList<>();
		for (AveragingProvider ap : averagingProviders) {
			List<Set<SpeciesVertex>> candidates = ap.findAveragingCandidates(reactionTimescales);
			for (Set<SpeciesVertex> candidate : candidates)
				if (!previousSubnetworkToAveragingProviderMap.containsKey(candidate))
					previousSubnetworkToAveragingProviderMap.put(candidate, ap);
//			allCandidatesMap.put(ap, candidates);
			allCandidates.addAll(candidates);
		}
		return allCandidates;
	}

	@Override
	public void reset() {
		previousSubnetworksToAverage = null;
		previousSubnetworkToAveragingProviderMap = null;
	}

	@Override
	public void resamplePreviouslyAveragedSubnetworks(double[] x,
			List<Set<SpeciesVertex>> subnetworksToAverage, List<Set<SpeciesVertex>> previousSubnetworksToAverage) {
		// Resample states that have been averaged before but are no longer averaged
		Set<SpeciesVertex> allAveragingSpecies = new HashSet<SpeciesVertex>();
		for (Set<SpeciesVertex> subnetwork : subnetworksToAverage)
			allAveragingSpecies.addAll(subnetwork);
		for (Set<SpeciesVertex> subnetworkSpecies : previousSubnetworksToAverage)
			if (!allAveragingSpecies.containsAll(subnetworkSpecies)) {
				resampleFromSteadyStateDistribution(x, subnetworkSpecies);
			}
	}

	@Override
	public void resampleFromSteadyStateDistribution(double[] x, Set<SpeciesVertex> subnetworkSpecies) {
		AveragingProvider ap = previousSubnetworkToAveragingProviderMap.get(subnetworkSpecies);
		ap.resampleFromSteadyStateDistribution(x, subnetworkSpecies);
	}

}
