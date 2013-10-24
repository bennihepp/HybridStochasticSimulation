package ch.ethz.khammash.hybridstochasticsimulation.averaging;

import com.google.common.base.Predicate;

public class MaximumSizeSubnetworkFilter implements Predicate<SubnetworkDescription> {

	private final int maxSize;

	public MaximumSizeSubnetworkFilter(final int maxSize) {
		this.maxSize = maxSize;
	}

	public boolean apply(SubnetworkDescription subnetwork) {
		return subnetwork.getSubnetworkSpecies().size() <= maxSize;
	}

}
