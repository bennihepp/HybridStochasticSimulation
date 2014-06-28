package ch.ethz.bhepp.hybridstochasticsimulation.outputs;

import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FiniteTrajectory;

public interface FiniteOutputMapper {

	OutputRecord mapToOutput(FiniteTrajectory tr);

}
