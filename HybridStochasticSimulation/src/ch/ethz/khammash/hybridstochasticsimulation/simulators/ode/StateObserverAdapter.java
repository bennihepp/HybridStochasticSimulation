package ch.ethz.khammash.hybridstochasticsimulation.simulators.ode;

import java.util.List;

import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;
import ch.ethz.khammash.ode.StateObserver;

public class StateObserverAdapter implements StateObserver {

	private PDMPModel model;
	private List<TrajectoryRecorder> trajectoryRecorders;
	private int lengthOfPrimaryState;
	private double[] completeState;

	public StateObserverAdapter(PDMPModel model, List<TrajectoryRecorder> trajectoryRecorders) {
		this.model = model;
		this.trajectoryRecorders = trajectoryRecorders;
	}

	@Override
	public void initialize(double t0, double[] x0, double t1) {
		double[] primaryState = model.computePrimaryState(t0, x0);
		lengthOfPrimaryState = primaryState.length;
		if (model.hasOptionalState()) {
			double[] optionalState = model.computeOptionalState(t0, x0);
			completeState = new double[lengthOfPrimaryState + optionalState.length];
			copyToFullState(primaryState, optionalState);
		} else
			completeState = primaryState;
		for (TrajectoryRecorder tr : trajectoryRecorders)
			tr.beginRecording(t0, completeState, t1);
	}

	@Override
	public void report(double t, double[] x) {
		double[] primaryState = model.computePrimaryState(t, x);
		if (model.hasOptionalState())
			copyToFullState(primaryState, model.computeOptionalState(t, x));
		else
			completeState = primaryState;
		for (TrajectoryRecorder tr : trajectoryRecorders)
			tr.record(t, completeState);
	}

	public void copyToFullState(double[] primaryState, double[] optionalState) {
		int i = 0;
		for (int s=0; s < lengthOfPrimaryState; s++)
			completeState[i++] = primaryState[s];
		for (int j=0; j < optionalState.length; j++)
			completeState[i++] = optionalState[j];
	}

}
