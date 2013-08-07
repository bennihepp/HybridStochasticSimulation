package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

public class DummyTrajectoryMapper implements FiniteTrajectoryMapper {

	@Override
	public FiniteTrajectory map(FiniteTrajectory tr) {
		return tr;
	}

}
