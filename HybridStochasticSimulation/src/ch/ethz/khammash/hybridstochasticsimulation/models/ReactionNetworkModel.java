package ch.ethz.khammash.hybridstochasticsimulation.models;

import ch.ethz.khammash.hybridstochasticsimulation.networks.ReactionNetwork;


public interface ReactionNetworkModel {

	int getNumberOfSpecies();

	int getNumberOfReactions();

	ReactionNetwork getNetwork();

}
