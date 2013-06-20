package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;


public class ArrayFiniteTrajectoryRecorder<T extends ReactionNetworkModel>
		extends ArrayFiniteTrajectory implements FiniteTrajectoryRecorder<T>, TrajectoryRecorder<T> {

	public ArrayFiniteTrajectoryRecorder(double[] tSeries) {
		super(tSeries);
	}

	public ArrayFiniteTrajectoryRecorder(double[] tSeries, double[][] xSeries) {
		super(tSeries, xSeries);
	}

	protected void initialize(double[] x0) {
		initialize(x0, x0.length);
	}

	protected void initialize(double[] x0, int numberOfStates) {
		xSeries = new double[numberOfStates][tSeries.length];
		index = 0;
	}

	protected void setState(int index, double[] x) {
		for (int s=0; s < xSeries.length; s++)
			xSeries[s][index] = x[s];
	}

	@Override
	public void setModel(T model) {
	}

	@Override
	public void setInitialState(double t0, double[] x0) {
		setInitialState(t0, x0, x0.length);
	}

	@Override
	public void setInitialState(double t0, double[] x0, int numOfStates) {
		initialize(x0, numOfStates);
		setState(index, x0);
	}

	@Override
	public void setFinalState(double t1, double[] x1) {
		setState(xSeries[0].length - 1, x1);
	}

	@Override
	public void reportState(double t, double[] x) {
		if (t <= tSeries[index]) {
			setState(index, x);
			if (t == tSeries[index])
				index++;
		} else if (t > tSeries[index]) {
			while (t > tSeries[index + 1]) {
				index++;
				for (int s=0; s < xSeries.length; s++)
					xSeries[s][index] = xSeries[s][index - 1];
			}
			setState(index + 1, x);
			index++;
		}
	}

//	@Override
//	public void handleReactionEvent(int reaction, double t, double[] newX) {
//		while (index < tSeries.length && tSeries[index] < t ) {
//			setState(index, previousX);
//			index++;
//		}
//		previousX = newX.clone();
//	}

}
