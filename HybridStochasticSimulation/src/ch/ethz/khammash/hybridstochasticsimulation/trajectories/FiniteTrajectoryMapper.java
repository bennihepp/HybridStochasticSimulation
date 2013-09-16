package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import java.io.Serializable;

public interface FiniteTrajectoryMapper extends Serializable {

	FiniteTrajectory map(FiniteTrajectory tr);

}
