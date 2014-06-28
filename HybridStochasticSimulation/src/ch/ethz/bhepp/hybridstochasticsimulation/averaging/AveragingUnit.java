package ch.ethz.bhepp.hybridstochasticsimulation.averaging;

import ch.ethz.bhepp.hybridstochasticsimulation.models.AdaptiveMSHRNModel;

import com.google.common.base.Predicate;

public interface AveragingUnit {

	void setSubnetworkEnumerator(SubnetworkEnumerator subnetworkEnumerator);

//	List<SubnetworkDescription> getSubnetworksToAverageAndResampleState(double t, double[] x, Predicate<SubnetworkDescription> filter);

	void reset();

	void updateAveraging(AdaptiveMSHRNModel model, double t, double[] x,
			SubnetworkDescription subnetwork);

	void updateSubnetworkState(AdaptiveMSHRNModel model, double t, double[] x,
			SubnetworkDescription subnetwork);

	void sampleSubnetworkState(double t, double[] x, SubnetworkDescription subnetwork);

	Predicate<SubnetworkDescription> getSubnetworkFilter();

}
