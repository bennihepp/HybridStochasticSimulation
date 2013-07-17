package ch.ethz.khammash.hybridstochasticsimulation.networks;

import java.util.List;
import java.util.Set;

import ch.ethz.khammash.hybridstochasticsimulation.sandbox.ReactionNetworkGraph;
import ch.ethz.khammash.hybridstochasticsimulation.sandbox.SpeciesVertex;

public interface AveragingProvider {

	void init(double theta, UnaryBinaryReactionNetwork network, ReactionNetworkGraph graph, Set<SpeciesVertex> importantSpecies);

	List<Set<SpeciesVertex>> getSubnetworksToAverageAndResampleState(double t, double[] x, double[] speciesTimescales);

}
