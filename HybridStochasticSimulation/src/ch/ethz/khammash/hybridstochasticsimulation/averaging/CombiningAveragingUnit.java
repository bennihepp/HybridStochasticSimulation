package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;

import com.google.common.base.Predicate;

public class CombiningAveragingUnit extends AbstractAveragingUnit {

	private List<AveragingUnit> averagingUnits;
	private List<Set<SpeciesVertex>> previousSubnetworksToAverage;
	private Map<Set<SpeciesVertex>, AveragingUnit> previousSubnetworkToAveragingUnitMap;

	public CombiningAveragingUnit() {
		averagingUnits = new LinkedList<>();
	}

	public void addAveragingUnit(AveragingUnit ap) {
		averagingUnits.add(ap);
	}

	public void removeAveragingUnit(AveragingUnit ap) {
		averagingUnits.remove(ap);
	}

	public void clearAveragingUnits() {
		averagingUnits.clear();
	}

	@Override
	public List<Set<SpeciesVertex>> getSubnetworksToAverageAndResampleState(double t, double[] x, Predicate<Set<SpeciesVertex>> filter) {
		List<Set<SpeciesVertex>> averagingCandidates = findAveragingCandidates(t, x, filter);
		List<Set<SpeciesVertex>> subnetworksToAverage = greedySelectSubnetworksToAverage(averagingCandidates);

		if (previousSubnetworksToAverage != null)
			resamplePreviouslyAveragedSubnetworks(t, x, subnetworksToAverage, previousSubnetworksToAverage);

		previousSubnetworksToAverage = subnetworksToAverage;
		return subnetworksToAverage;
	}

	@Override
	public List<Set<SpeciesVertex>> findAveragingCandidates(double t, double[] x, Predicate<Set<SpeciesVertex>> filter) {
		previousSubnetworkToAveragingUnitMap = new HashMap<>();
		List<Set<SpeciesVertex>> allCandidates = new ArrayList<>();
		for (AveragingUnit au : averagingUnits) {
			List<Set<SpeciesVertex>> candidates = au.findAveragingCandidates(t, x, filter);
			for (Set<SpeciesVertex> candidate : candidates)
				if (!previousSubnetworkToAveragingUnitMap.containsKey(candidate))
					previousSubnetworkToAveragingUnitMap.put(candidate, au);
			allCandidates.addAll(candidates);
		}
		return allCandidates;
	}

	@Override
	public void reset() {
		previousSubnetworksToAverage = null;
		previousSubnetworkToAveragingUnitMap = null;
	}

	@Override
	public void resamplePreviouslyAveragedSubnetworks(double t, double[] x,
			List<Set<SpeciesVertex>> subnetworksToAverage, List<Set<SpeciesVertex>> previousSubnetworksToAverage) {
		// Resample states that have been averaged before but are no longer averaged
		Set<SpeciesVertex> allAveragingSpecies = new HashSet<SpeciesVertex>();
		for (Set<SpeciesVertex> subnetwork : subnetworksToAverage)
			allAveragingSpecies.addAll(subnetwork);
		for (Set<SpeciesVertex> subnetworkSpecies : previousSubnetworksToAverage)
			if (!allAveragingSpecies.containsAll(subnetworkSpecies)) {
				resampleFromSteadyStateDistribution(t, x, subnetworkSpecies);
			}
	}

	@Override
	public void resampleFromSteadyStateDistribution(double t, double[] x, Set<SpeciesVertex> subnetworkSpecies) {
		AveragingUnit ap = previousSubnetworkToAveragingUnitMap.get(subnetworkSpecies);
		ap.resampleFromSteadyStateDistribution(t, x, subnetworkSpecies);
	}

}
