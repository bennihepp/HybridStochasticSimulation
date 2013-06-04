package ch.ethz.khammash.hybridstochasticsimulation.simulators;

import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;

public class PDMPFixedModelTrajectory implements PDMPStepHandler, ReactionEventHandler {

	double[] tSeries;
	double[][] xSeries;
	int index;
	private boolean initialized;

	public PDMPFixedModelTrajectory(double[] tSeries) {
		this.tSeries = tSeries;
		initialized = false;
	}

	private void initialize(double t0, double[] x0) {
		if (initialized == false) {
			xSeries = new double[x0.length - 2][tSeries.length];
			for (int i=0; i < xSeries.length; i++)
				xSeries[i][0] = x0[i];
			initialized = true;
			index = 0;
		}
	}

	public void reset() {
		initialized = false;
	}

	@Override
	public void init(double t0, double[] x0, double t) {
		initialize(t0, x0);
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
			for (int i=0; i < xSeries.length; i++)
				xSeries[i][index] = x[i];
			index++;
		}
	}

	@Override
	public void setInitialState(double t0, double[] x0) {
		initialize(t0, x0);
	}

	@Override
	public void setFinalState(double t1, double[] x1) {
	}

	@Override
	public void handleReactionEvent(int reaction, double t, double[] newX) {
		while (index < tSeries.length && t >= tSeries[index]) {
			for (int i=0; i < xSeries.length; i++)
				xSeries[i][index] = newX[i];
			index++;
		}
	}

	@Override
	public void setPDMPModel(PDMPModel model) {
	}

}
