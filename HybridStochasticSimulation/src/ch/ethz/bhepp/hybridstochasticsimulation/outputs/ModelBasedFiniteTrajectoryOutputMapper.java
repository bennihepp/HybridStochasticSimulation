package ch.ethz.bhepp.hybridstochasticsimulation.outputs;

import java.util.List;

import javax.inject.Inject;

import ch.ethz.bhepp.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FiniteTrajectory;

public class ModelBasedFiniteTrajectoryOutputMapper<T extends ReactionNetworkModel> extends ModelBasedOutputMapper<T> {

	@Inject
	public ModelBasedFiniteTrajectoryOutputMapper(List<String> speciesLabels) {
		super(speciesLabels);
	}

	@Override
	public OutputRecord mapToOutput(FiniteTrajectory tr) {
		return new FiniteTrajectoryOutput(tr, getLabels());
	}

}
