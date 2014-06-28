package ch.ethz.bhepp.hybridstochasticsimulation.trajectories;

public class DummyTrajectoryMapper implements FiniteTrajectoryMapper {

	private static final long serialVersionUID = 1L;

	@Override
	public FiniteTrajectory map(FiniteTrajectory tr) {
		return tr;
	}

}
