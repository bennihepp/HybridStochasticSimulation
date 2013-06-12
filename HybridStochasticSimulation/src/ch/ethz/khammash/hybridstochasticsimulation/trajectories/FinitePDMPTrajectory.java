package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;

public class FinitePDMPTrajectory<T extends ReactionNetworkModel>
		extends FiniteStochasticTrajectory<T> implements ContinuousTrajectoryRecorder<T> {

	public FinitePDMPTrajectory(double[] tSeries) {
		super(tSeries);
	}

	@Override
	public void init(double t0, double[] x0, double t) {
	}

	@Override
	public void handleStep(StepInterpolator interpolator, boolean isLast)
			throws MaxCountExceededException {
		while (index < tSeries.length && interpolator.getCurrentTime() >= tSeries[index]) {
			interpolator.setInterpolatedTime(tSeries[index]);
			double[] x = interpolator.getInterpolatedState();
			setState(index, x);
			index++;
		}
	}

}
