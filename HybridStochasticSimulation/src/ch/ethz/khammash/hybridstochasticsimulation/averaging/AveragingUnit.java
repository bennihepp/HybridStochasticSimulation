package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import java.util.List;
import java.util.Set;

import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;

public interface AveragingUnit {

	List<Set<SpeciesVertex>> getSubnetworksToAverageAndResampleState(double t, double[] x, double[] reactionTimescales);

	List<Set<SpeciesVertex>> findAveragingCandidates(double t, double[] x, double[] reactionTimescales);

	void reset();

	void resampleFromSteadyStateDistribution(double t, double[] x, Set<SpeciesVertex> subnetworkSpecies);

}