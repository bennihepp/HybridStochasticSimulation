package ch.ethz.khammash.hybridstochasticsimulation.outputs;

import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;

public interface FiniteOutputMapper {

	OutputRecord mapToOutput(FiniteTrajectory tr);

}
