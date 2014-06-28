package ch.ethz.bhepp.hybridstochasticsimulation.averaging;

import javax.inject.Inject;

import ch.ethz.bhepp.hybridstochasticsimulation.models.AdaptiveMSHRNModel;

import com.google.common.base.Predicate;

public class DummyAveragingUnit implements AveragingUnit {

	@Inject
	public DummyAveragingUnit() {}

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

	@Override
	public Predicate<SubnetworkDescription> getSubnetworkFilter() {
		return new Predicate<SubnetworkDescription>() {

			@Override
			public boolean apply(SubnetworkDescription subnetworkDescr) {
				return false;
			}

		};
	}

}
