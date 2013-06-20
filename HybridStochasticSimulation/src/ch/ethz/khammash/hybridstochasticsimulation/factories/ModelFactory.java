package ch.ethz.khammash.hybridstochasticsimulation.factories;

import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;

public interface ModelFactory<T extends ReactionNetworkModel> {

	T createModel();

}
