package ch.ethz.khammash.hybridstochasticsimulation.models;

import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPStepHandler;

public class FinitePDMPModelTrajectory extends FiniteStochasticModelTrajectory implements PDMPStepHandler {

	protected boolean initialized;

	public FinitePDMPModelTrajectory(double[] tSeries) {
		super(tSeries);
		initialized = false;
	}

	@Override
	public void reset() {
		initialized = false;
	}

	@Override
	public void setPDMPModel(PDMPModel model) {
	}

	@Override
	protected void initialize(double[] x0, int numberOfStates) {
		if (initialized == false) {
			super.initialize(x0, numberOfStates);
			initialized = true;
		}
	}

	@Override
	public void init(double t0, double[] x0, double t) {
		initialize(x0, x0.length - 2);
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
