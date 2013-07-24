package ch.ethz.khammash.hybridstochasticsimulation.mpj;

import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;

public class SimulationResultContainer {

	private int mpiSource;
	private FiniteTrajectory tr;

	public SimulationResultContainer(int mpiSource, FiniteTrajectory tr) {
		this.mpiSource = mpiSource;
		this.tr = tr;
	}

	public int getMpiSource() {
		return mpiSource;
	}

	public FiniteTrajectory getFiniteTrajectory() {
		return tr;
	}

}
