package ch.ethz.khammash.hybridstochasticsimulation.batch;

import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;

public interface ClientCallbackInterface {

	void callback(FiniteTrajectory tr);

}
