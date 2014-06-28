package ch.ethz.bhepp.hybridstochasticsimulation.models;

import ch.ethz.bhepp.hybridstochasticsimulation.networks.MassActionReactionNetwork;


public interface ReactionNetworkModel {

	int getNumberOfSpecies();

	int getNumberOfReactions();

	MassActionReactionNetwork getNetwork();

}
