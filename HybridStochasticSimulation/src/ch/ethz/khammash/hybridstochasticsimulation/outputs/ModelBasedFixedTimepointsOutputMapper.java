package ch.ethz.khammash.hybridstochasticsimulation.outputs;

import java.util.List;

import javax.inject.Inject;
import javax.inject.Named;

import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;

public class ModelBasedFixedTimepointsOutputMapper<T extends ReactionNetworkModel> extends ModelBasedOutputMapper<T> {

	private double time;

	@Inject
	public ModelBasedFixedTimepointsOutputMapper(@Named("outputTimepoints") double time, @Named("speciesLabels") List<String> speciesLabels) {
		super(speciesLabels);
		this.time = time;
	}

	@Override
	public OutputRecord mapToOutput(FiniteTrajectory tr) {
		return new FixedTimepointOutput(time, tr, getLabels());
	}

}
