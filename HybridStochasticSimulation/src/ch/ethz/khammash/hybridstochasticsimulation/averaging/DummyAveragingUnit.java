package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import java.util.Collections;
import java.util.List;
import java.util.Set;

import javax.inject.Inject;

import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;

public class DummyAveragingUnit implements AveragingUnit {

	@Inject
	public DummyAveragingUnit() {}

	@Override
	public List<Set<SpeciesVertex>> getSubnetworksToAverageAndResampleState(double t, double[] x, double[] reactionTimescales) {
		return Collections.<Set<SpeciesVertex>>emptyList();
	}

	@Override
	public List<Set<SpeciesVertex>> findAveragingCandidates(double t, double[] x, double[] reactionTimescales) {
		return Collections.<Set<SpeciesVertex>>emptyList();
	}

	@Override
	public void reset() {
	}

	@Override
	public void resampleFromSteadyStateDistribution(double t, double[] x, Set<SpeciesVertex> subnetworkSpecies) {
	}

}
