package ch.ethz.bhepp.hybridstochasticsimulation.grid.simon;

import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FiniteTrajectory;

public interface ControllerInterface {

	boolean moreSimulations();

	void addSimulationResult(FiniteTrajectory tr);

	void register(Object object);

	void unregister(Object object);

}
