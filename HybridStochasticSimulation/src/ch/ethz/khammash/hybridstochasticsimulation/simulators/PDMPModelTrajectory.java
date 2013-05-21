package ch.ethz.khammash.hybridstochasticsimulation.simulators;

import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.sampling.StepInterpolator;


public class PDMPModelTrajectory extends StochasticModelTrajectory implements ModelTrajectory, PDMPStepHandler, ReactionHandler {

	ContinuousOutputModel com;
	protected boolean initialized;

	public PDMPModelTrajectory() {
		com = new ContinuousOutputModel();
		initialized = false;
	}

	public void reset() {
		initialized = false;
	}

	@Override
	public void init(double t0, double[] y0, double t) {
		if (initialized == false) {
			initialized = true;
			com.init(t0, y0, t);
		}
	}

	@Override
	public double getInitialtime() {
		return com.getInitialTime();
	}

	@Override
	public double getFinalTime() {
		return com.getFinalTime();
	}

	public void setInterpolatedTime(double t) {
		com.setInterpolatedTime(t);
	}

	public double[] getInterpolatedState() {
		double[] x = com.getInterpolatedState();
		double[] reducedX = new double[x.length - 2];
		for (int i = 0; i < reducedX.length; i++)
			reducedX[i] = x[i];
		return reducedX;
	}

	@Override
	public double[] getInterpolatedState(double t) {
		setInterpolatedTime(t);
		return getInterpolatedState();
	}

	@Override
	public void handleStep(StepInterpolator interpolator, boolean isLast) throws MaxCountExceededException {
		com.handleStep(interpolator, isLast);
	}

}
