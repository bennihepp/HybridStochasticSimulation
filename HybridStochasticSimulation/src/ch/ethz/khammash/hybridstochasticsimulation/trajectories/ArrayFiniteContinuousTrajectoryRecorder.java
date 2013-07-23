package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

public class ArrayFiniteContinuousTrajectoryRecorder extends ArrayFiniteTrajectoryRecorder implements ContinuousTrajectoryRecorder {

	private static final long serialVersionUID = 6341601501856273860L;

	private int index;

	public ArrayFiniteContinuousTrajectoryRecorder(int numOfTimePoints) {
		super(numOfTimePoints);
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
