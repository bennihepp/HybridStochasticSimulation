package ch.ethz.khammash.hybridstochasticsimulation.grid.simon;

import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;

public interface ControllerInterface {

	boolean moreSimulations();

	void addSimulationResult(FiniteTrajectory tr);

	void register(Object object);

	void unregister(Object object);

}
