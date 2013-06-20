package ch.ethz.khammash.hybridstochasticsimulation.simulators.lsodar;

import java.util.List;

import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;
import ch.ethz.khammash.nativeode.StateObserver;

public class LsodarStateObserverAdapter<T extends PDMPModel> implements StateObserver {

	private List<FiniteTrajectoryRecorder<T>> trajectoryRecorders;

	public LsodarStateObserverAdapter(List<FiniteTrajectoryRecorder<T>> trajectoryRecorders) {
		this.trajectoryRecorders = trajectoryRecorders;
	}

	@Override
	public void initialize(double t0, double[] x0, double t1) {
		for (FiniteTrajectoryRecorder<T> tr : trajectoryRecorders)
			tr.setInitialState(t0, x0);
	}

	@Override
	public void report(double t, double[] x) {
		for (FiniteTrajectoryRecorder<T> tr : trajectoryRecorders)
			tr.reportState(t, x);
	}

}
