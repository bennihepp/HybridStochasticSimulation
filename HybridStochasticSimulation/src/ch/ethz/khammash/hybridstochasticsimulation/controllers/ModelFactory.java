package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;

public interface ModelFactory<T extends ReactionNetworkModel> {

	public T createModel();

}
