package ch.ethz.khammash.hybridstochasticsimulation.hazelcast;

import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;

public interface ResultHandler {

	void handle(FiniteTrajectory tr);

}
