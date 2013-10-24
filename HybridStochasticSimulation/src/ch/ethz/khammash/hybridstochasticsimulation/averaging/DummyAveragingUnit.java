package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import java.util.Collections;
import java.util.List;

import javax.inject.Inject;

import ch.ethz.khammash.hybridstochasticsimulation.models.AdaptiveMSHRNModel;

import com.google.common.base.Predicate;

public class DummyAveragingUnit implements AveragingUnit {

	@Inject
	public DummyAveragingUnit() {}

	@Override
	public List<SubnetworkDescription> getSubnetworksToAverageAndResampleState(double t, double[] x, Predicate<SubnetworkDescription> filter) {
		return Collections.emptyList();
	}

	@Override
	public void reset() {
	}

	@Override
	public void setSubnetworkEnumerator(SubnetworkEnumerator subnetworksEnumerator) {
	}

	@Override
	public void updateAveraging(AdaptiveMSHRNModel model, double t, double[] x,
			SubnetworkDescription subnetwork) {
		throw new UnsupportedOperationException();
	}

	@Override
	public void sampleSubnetworkState(double t, double[] x,
			SubnetworkDescription subnetwork) {
		throw new UnsupportedOperationException();
	}

	@Override
	public void updateSubnetworkState(AdaptiveMSHRNModel model, double t,
			double[] x, SubnetworkDescription subnetwork) {
		throw new UnsupportedOperationException();
	}

}
