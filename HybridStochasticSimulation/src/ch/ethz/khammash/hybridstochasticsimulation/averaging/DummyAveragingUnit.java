package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import java.util.Collections;
import java.util.List;
import java.util.Set;

import javax.inject.Inject;

import ch.ethz.khammash.hybridstochasticsimulation.graphs.SpeciesVertex;

import com.google.common.base.Predicate;

public class DummyAveragingUnit implements AveragingUnit {

	@Inject
	public DummyAveragingUnit() {}

	@Override
	public List<Set<SpeciesVertex>> getSubnetworksToAverageAndResampleState(double t, double[] x, Predicate<Set<SpeciesVertex>> filter) {
		return Collections.<Set<SpeciesVertex>>emptyList();
	}

	@Override
	public void reset() {
	}

	@Override
	public void setSubnetworkEnumerator(SubnetworkEnumerator subnetworksEnumerator) {
	}

}
