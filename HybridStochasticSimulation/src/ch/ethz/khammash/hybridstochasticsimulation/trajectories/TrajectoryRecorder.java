package ch.ethz.khammash.hybridstochasticsimulation.trajectories;


public interface TrajectoryRecorder extends Trajectory {

	void beginRecording(double t0, double[] x0, double t1);

	void endRecording(double[] x1);

//	void setInitialState(double t0, double[] x0);

//	void setInitialState(double t0, double[] x0, int numOfStates);

//	void setFinalState(double t1, double[] x1);

	void record(double t, double[] x);

}
