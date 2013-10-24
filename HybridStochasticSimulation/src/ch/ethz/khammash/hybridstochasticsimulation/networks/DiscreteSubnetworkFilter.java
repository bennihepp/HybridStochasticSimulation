package ch.ethz.khammash.hybridstochasticsimulation.networks;

import ch.ethz.khammash.hybridstochasticsimulation.averaging.SubnetworkDescription;
import ch.ethz.khammash.hybridstochasticsimulation.networks.MSHybridReactionNetwork.ReactionType;

import com.google.common.base.Predicate;

public class DiscreteSubnetworkFilter implements Predicate<SubnetworkDescription> {

	private MSHybridReactionNetwork hrn;

	public DiscreteSubnetworkFilter(MSHybridReactionNetwork hrn) {
		this.hrn = hrn;
	}

	@Override
	public boolean apply(SubnetworkDescription subnetwork) {
//		for (int s : subnetwork.getSubnetworkSpeciesIndices()) {
//			if (hrn.getSpeciesType(s) != SpeciesType.DISCRETE)
//				return false;
//		}
		for (int r : subnetwork.getSubnetworkReactions()) {
			if (hrn.getReactionType(r) != ReactionType.DISCRETE)
				return false;
		}
		return true;
	}

}
