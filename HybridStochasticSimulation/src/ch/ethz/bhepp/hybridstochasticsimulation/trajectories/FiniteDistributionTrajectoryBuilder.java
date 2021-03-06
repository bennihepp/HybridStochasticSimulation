package ch.ethz.bhepp.hybridstochasticsimulation.trajectories;

public interface FiniteDistributionTrajectoryBuilder {

	void addTrajectory(FiniteTrajectory tr);

	FiniteDistributionTrajectory getDistributionTrajectory();

	int getNumberOfAddedTrajectories();

}
