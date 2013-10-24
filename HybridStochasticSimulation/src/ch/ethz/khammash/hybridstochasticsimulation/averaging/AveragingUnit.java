package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import java.util.List;

import ch.ethz.khammash.hybridstochasticsimulation.models.AdaptiveMSHRNModel;

import com.google.common.base.Predicate;

public interface AveragingUnit {

	void setSubnetworkEnumerator(SubnetworkEnumerator subnetworkEnumerator);

	List<SubnetworkDescription> getSubnetworksToAverageAndResampleState(double t, double[] x, Predicate<SubnetworkDescription> filter);

	void reset();

	void updateAveraging(AdaptiveMSHRNModel model, double t, double[] x,
			SubnetworkDescription subnetwork);

	void updateSubnetworkState(AdaptiveMSHRNModel model, double t, double[] x,
			SubnetworkDescription subnetwork);

	void sampleSubnetworkState(double t, double[] x, SubnetworkDescription subnetwork);

}
