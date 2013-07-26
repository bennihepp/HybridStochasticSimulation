package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import java.util.Set;

import ch.ethz.khammash.hybridstochasticsimulation.graphs.ReactionNetworkGraph;
import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;

import com.google.common.base.Predicate;

public class MaximumSizeSubnetworkEnumerator extends FilteredSubnetworksEnumerator {

	public MaximumSizeSubnetworkEnumerator(ReactionNetworkGraph graph, final int maxSize) {
		super(graph, new Predicate<Set<SpeciesVertex>>() {

			public boolean apply(Set<SpeciesVertex> input) {
				return input.size() <= maxSize;
			}

		});
	}

}
