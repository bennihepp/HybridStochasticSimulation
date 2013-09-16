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

public class CombiningAveragingUnit implements AveragingUnit {

	private List<ModularAveragingUnit> averagingUnits;
	private List<Set<SpeciesVertex>> lastSubnetworksToAverage;
	private Map<Set<SpeciesVertex>, ModularAveragingUnit> lastSubnetworkToAveragingUnitMap;
	private SubnetworkEnumerator subnetworkEnumerator;

	public CombiningAveragingUnit() {
		averagingUnits = new LinkedList<>();
	}

	public void addAveragingUnit(ModularAveragingUnit au) {
		au.setSubnetworkEnumerator(subnetworkEnumerator);
		averagingUnits.add(au);
	}

	public void removeAveragingUnit(ModularAveragingUnit ap) {
		averagingUnits.remove(ap);
	}

	public void clearAveragingUnits() {
		averagingUnits.clear();
	}

	@Override
	public void setSubnetworkEnumerator(SubnetworkEnumerator subnetworkEnumerator) {
		this.subnetworkEnumerator = subnetworkEnumerator;
		for (ModularAveragingUnit au : averagingUnits)
			au.setSubnetworkEnumerator(subnetworkEnumerator);
		
	}

	@Override
	public List<Set<SpeciesVertex>> getSubnetworksToAverageAndResampleState(double t, double[] x, Predicate<Set<SpeciesVertex>> filter) {
		Map<Set<SpeciesVertex>, ModularAveragingUnit> subnetworkToAveragingUnitMap = findAveragingCandidates(t, x, filter);
		List<Set<SpeciesVertex>> averagingCandidates = new ArrayList<>(subnetworkToAveragingUnitMap.keySet());
		List<Set<SpeciesVertex>> subnetworksToAverage = AbstractAveragingUnit.greedySelectSubnetworksToAverage(averagingCandidates);

		if (lastSubnetworksToAverage != null && lastSubnetworksToAverage.size() > 0)
			resampleLastAveragedSubnetworks(t, x, subnetworksToAverage);

		lastSubnetworksToAverage = subnetworksToAverage;
		lastSubnetworkToAveragingUnitMap = subnetworkToAveragingUnitMap;

		return subnetworksToAverage;
	}

	private Map<Set<SpeciesVertex>, ModularAveragingUnit> findAveragingCandidates(double t, double[] x, Predicate<Set<SpeciesVertex>> filter) {
		Map<Set<SpeciesVertex>, ModularAveragingUnit> subnetworkToAveragingUnitMap = new HashMap<>();
//		List<Set<SpeciesVertex>> allCandidates = new ArrayList<>();
		for (ModularAveragingUnit au : averagingUnits) {
			List<Set<SpeciesVertex>> candidates = au.findAveragingCandidates(t, x, filter);
			for (Set<SpeciesVertex> candidate : candidates)
				if (!subnetworkToAveragingUnitMap.containsKey(candidate)) {
					subnetworkToAveragingUnitMap.put(candidate, au);
//					allCandidates.add(candidate);
				}
		}
		return subnetworkToAveragingUnitMap;
//		return allCandidates;
	}

	@Override
	public void reset() {
		lastSubnetworksToAverage = null;
		lastSubnetworkToAveragingUnitMap = null;
	}

	private void resampleLastAveragedSubnetworks(double t, double[] x, List<Set<SpeciesVertex>> subnetworksToAverage) {
		// Resample states that have been averaged before but are no longer averaged
		Set<SpeciesVertex> allAveragingSpecies = new HashSet<SpeciesVertex>();
		for (Set<SpeciesVertex> subnetwork : subnetworksToAverage)
			allAveragingSpecies.addAll(subnetwork);
		for (Set<SpeciesVertex> subnetworkSpecies : lastSubnetworksToAverage)
			if (!allAveragingSpecies.containsAll(subnetworkSpecies)) {
				resampleFromSteadyStateDistribution(t, x, subnetworkSpecies);
			}
	}

	private void resampleFromSteadyStateDistribution(double t, double[] x, Set<SpeciesVertex> subnetworkSpecies) {
		ModularAveragingUnit ap = lastSubnetworkToAveragingUnitMap.get(subnetworkSpecies);
		ap.resampleFromStationaryDistribution(t, x, subnetworkSpecies);
	}

}
