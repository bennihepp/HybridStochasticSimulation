package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import org.apache.commons.math3.ode.sampling.StepHandler;

import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;

public interface ContinuousTrajectoryRecorder<T extends ReactionNetworkModel> extends TrajectoryRecorder<T>, StepHandler {

}
