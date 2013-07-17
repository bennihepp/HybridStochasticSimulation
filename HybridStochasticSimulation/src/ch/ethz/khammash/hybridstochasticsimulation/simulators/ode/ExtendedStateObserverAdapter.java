package ch.ethz.khammash.hybridstochasticsimulation.simulators.ode;

import java.util.List;

import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.SimulationInformation;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;

public class ExtendedStateObserverAdapter extends StateObserverAdapter {

	private SimulationInformation simInfo;
	private List<TrajectoryRecorder> simulationInformationTrajectoryRecorders;

	public ExtendedStateObserverAdapter(SimulationInformation simInfo, PDMPModel model,
			List<TrajectoryRecorder> trajectoryRecorders, List<TrajectoryRecorder> optionalTrajectoryRecorders,
			List<TrajectoryRecorder> simulationInformationTrajectoryRecorders) {
		super(model, trajectoryRecorders, optionalTrajectoryRecorders);
		this.simInfo = simInfo;
		this.simulationInformationTrajectoryRecorders = simulationInformationTrajectoryRecorders;
	}

	@Override
	public void initialize(double t0, double[] x0, double t1) {
		super.initialize(t0, x0, t1);
		double[] simInfoState = simInfo.computeInformationState();
		for (TrajectoryRecorder tr : simulationInformationTrajectoryRecorders)
			tr.beginRecording(t0, simInfoState, t1);
	}

	@Override
	public void report(double t, double[] x) {
		super.report(t, x);
		double[] simInfoState = simInfo.computeInformationState();
		for (TrajectoryRecorder tr : simulationInformationTrajectoryRecorders)
			tr.record(t, simInfoState);
	}

//	@Override
//	protected double[] computeOptionalState(double t, double[] x) {
//		double[] optionalState = super.computeOptionalState(t, x);
//		double[] simInfoState = simInfo.getInformationState();
//		double[] newOptionalState = new double[optionalState.length + simInfoState.length];
//		System.arraycopy(optionalState, 0, newOptionalState, 0, optionalState.length);
//		System.arraycopy(simInfoState, 0, newOptionalState, optionalState.length, simInfoState.length);
////		int i = optionalState.length;
////		newOptionalState[i++] = simInfo.getIntegrationCount();
////		newOptionalState[i++] = simInfo.getReactionCount();
////		newOptionalState[i++] = simInfo.isIntegrating() ? 1.0 : -1.0;
//		return newOptionalState;
//	}

}
