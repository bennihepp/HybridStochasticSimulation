package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;


public class FiniteStochasticTrajectory<T extends ReactionNetworkModel> extends FiniteTrajectory implements TrajectoryRecorder<T> {

	public FiniteStochasticTrajectory(double[] tSeries) {
		super(tSeries);
	}

	@Override
	public void setModel(T model) {
	}

	@Override
	public void setInitialState(double t0, double[] x0) {
		initialize(x0, x0.length);
		setState(index, x0);
	}

	@Override
	public void setFinalState(double t1, double[] x1) {
		setState(xSeries[0].length - 1, x1);
	}

	@Override
	public void handleReactionEvent(int reaction, double t, double[] newX) {
		while (index < tSeries.length && t >= tSeries[index]) {
			setState(index, newX);
			index++;
		}
	}

}
