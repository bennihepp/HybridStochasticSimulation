package ch.ethz.bhepp.hybridstochasticsimulation.outputs;

import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.Trajectory;

public interface OutputMapper {

	OutputRecord mapToOutput(Trajectory tr);

}
