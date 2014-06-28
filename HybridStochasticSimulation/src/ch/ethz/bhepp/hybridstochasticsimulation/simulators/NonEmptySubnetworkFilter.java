package ch.ethz.bhepp.hybridstochasticsimulation.simulators;

import ch.ethz.bhepp.hybridstochasticsimulation.averaging.SubnetworkDescription;

import com.google.common.base.Predicate;

public class NonEmptySubnetworkFilter implements Predicate<SubnetworkDescription> {

	@Override
	public boolean apply(SubnetworkDescription subnetwork) {
		return !subnetwork.getSubnetworkSpecies().isEmpty()
				&& !subnetwork.getSubnetworkReactions().isEmpty()
				&& !subnetwork.getSurroundingReactions().isEmpty();
	}

}
