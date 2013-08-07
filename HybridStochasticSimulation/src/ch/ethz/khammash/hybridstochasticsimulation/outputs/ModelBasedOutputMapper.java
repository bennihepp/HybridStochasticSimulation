package ch.ethz.khammash.hybridstochasticsimulation.outputs;

import java.util.Collections;
import java.util.List;

import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;

public abstract class ModelBasedOutputMapper<T extends ReactionNetworkModel> implements FiniteOutputMapper {

	private List<String> labels;

	public ModelBasedOutputMapper(List<String> labels) {
		this.labels = labels;
	}

	public List<String> getLabels() {
		return Collections.unmodifiableList(labels);
	}

}
