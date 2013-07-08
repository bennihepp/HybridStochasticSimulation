package ch.ethz.khammash.hybridstochasticsimulation.trajectories;


public class ArrayFiniteTrajectoryRecorder extends ArrayFiniteTrajectory implements FiniteTrajectoryRecorder, TrajectoryRecorder {

	private int numOfTimePoints;

	public ArrayFiniteTrajectoryRecorder(int numOfTimePoints) {
		super(new double[numOfTimePoints]);
		this.numOfTimePoints = numOfTimePoints;
	}

	protected void initialize(double t0, double[] x0, double t1) {
		tSeries = computeTimeSeries(numOfTimePoints, t0, t1);
		xSeries = new double[x0.length][tSeries.length];
		index = 0;
	}

	public static double[] computeTimeSeries(int numOfTimePoints, double t0, double t1) {
		double[] tSeries = new double[numOfTimePoints];
		for (int i = 0; i < numOfTimePoints; i++) {
			tSeries[i] = t0 + i * (t1 - t0) / (double) (numOfTimePoints - 1);
		}
		return tSeries;
	}

	protected void setState(int index, double[] x) {
		for (int s=0; s < xSeries.length; s++)
			xSeries[s][index] = x[s];
	}

	@Override
	public void beginRecording(double t0, double[] x0, double t1) {
		initialize(t0, x0, t1);
		setState(index, x0);
	}

	@Override
	public void endRecording(double[] x1) {
		setState(xSeries[0].length - 1, x1);
	}

	@Override
	public void record(double t, double[] x) {
		if (index >= tSeries.length)
			return;
		if (t <= tSeries[index]) {
			setState(index, x);
			if (t == tSeries[index])
				index++;
		} else {
			while (index + 1 < tSeries.length && t > tSeries[index + 1]) {
				index++;
				for (int s=0; s < xSeries.length; s++)
					xSeries[s][index] = xSeries[s][index - 1];
			}
			index++;
			if (index < tSeries.length)
				setState(index, x);
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
