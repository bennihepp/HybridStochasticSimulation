package ch.ethz.bhepp.hybridstochasticsimulation.simulators.ode;

import java.util.List;

import ch.ethz.bhepp.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.TrajectoryRecorder;
import ch.ethz.bhepp.ode.StateObserver;

public class StateObserverAdapter implements StateObserver {

	private PDMPModel model;
	private List<TrajectoryRecorder> trajectoryRecorders;
	private List<TrajectoryRecorder> optionalTrajectoryRecorders;

	public StateObserverAdapter(PDMPModel model, List<TrajectoryRecorder> trajectoryRecorders, List<TrajectoryRecorder> optionalTrajectoryRecorders) {
		this.model = model;
		this.trajectoryRecorders = trajectoryRecorders;
		this.optionalTrajectoryRecorders = optionalTrajectoryRecorders;
	}

	@Override
	public void initialize(double t0, double[] x0, double t1) {
		double[] primaryState = computePrimaryState(t0, x0);
		for (TrajectoryRecorder tr : trajectoryRecorders)
			tr.beginRecording(t0, primaryState, t1);
		if (model.hasOptionalState()) {
			double[] optionalState = computeOptionalState(t0, x0);
			for (TrajectoryRecorder tr : optionalTrajectoryRecorders)
				tr.beginRecording(t0, optionalState, t1);
		}
	}

	@Override
	public void report(double t, double[] x) {
		double[] primaryState = computePrimaryState(t, x);
		for (TrajectoryRecorder tr : trajectoryRecorders)
			tr.record(t, primaryState);
		if (model.hasOptionalState()) {
			double[] optionalState = computeOptionalState(t, x);
			for (TrajectoryRecorder tr : optionalTrajectoryRecorders)
				tr.record(t, optionalState);
		}
	}

	protected double[] computePrimaryState(double t, double[] x) {
		return model.computePrimaryState(t, x);
	}

	protected double[] computeOptionalState(double t, double[] x) {
		return model.computeOptionalState(t, x);
	}

}
