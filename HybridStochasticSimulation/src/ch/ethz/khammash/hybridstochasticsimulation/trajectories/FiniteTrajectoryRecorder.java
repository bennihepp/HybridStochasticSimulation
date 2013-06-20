package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;

public interface FiniteTrajectoryRecorder<T extends ReactionNetworkModel> extends FiniteTrajectory, TrajectoryRecorder<T> {

}
