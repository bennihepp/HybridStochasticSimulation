package ch.ethz.khammash.hybridstochasticsimulation.models;

import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPStepHandler;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.ReactionEventHandler;

public class PDMPFixedModelTrajectory implements PDMPStepHandler, ReactionEventHandler {

	protected double[] tSeries;
	protected double[][] xSeries;
	protected int index;
	protected boolean initialized;

	public PDMPFixedModelTrajectory(double[] tSeries) {
		this.tSeries = tSeries;
		initialized = false;
	}

	@Override
	public void reset() {
		initialized = false;
	}

	@Override
	public void setPDMPModel(PDMPModel model) {
	}

	protected void initialize(double[] x0, int numberOfStates) {
		if (initialized == false) {
			xSeries = new double[numberOfStates][tSeries.length];
			setState(0, x0);
			initialized = true;
			index = 0;
		}
	}

	protected void setState(int index, double[] x) {
		for (int i=0; i < xSeries.length; i++)
			xSeries[i][index] = x[i];
	}

	@Override
	public void init(double t0, double[] x0, double t) {
		initialize(x0, x0.length - 2);
	}

	public double[] gettSeries() {
		return tSeries;
	}

	public double[][] getxSeries() {
		return xSeries;
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

	@Override
	public void setInitialState(double t0, double[] x0) {
		initialize(x0, x0.length);
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
