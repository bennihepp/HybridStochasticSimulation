package ch.ethz.khammash.hybridstochasticsimulation.outputs;

import ch.ethz.khammash.hybridstochasticsimulation.trajectories.Trajectory;

public interface OutputMapper {

	OutputRecord mapToOutput(Trajectory tr);

}
