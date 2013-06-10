package ch.ethz.khammash.hybridstochasticsimulation.models;

import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

public class FiniteDeterministicModelTrajectory extends FiniteModelTrajectory implements StepHandler {

	public FiniteDeterministicModelTrajectory(double[] tSeries) {
		super(tSeries);
	}

	@Override
	public void init(double t0, double[] x0, double t) {
		initialize(x0, x0.length);
		setState(0, x0);
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
