package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.ContinuousTrajectoryRecorder;

public interface FiniteContinuousTrajectoryRecorder<T extends ReactionNetworkModel>
	extends FiniteTrajectoryRecorder<T>, ContinuousTrajectoryRecorder<T> {

}
