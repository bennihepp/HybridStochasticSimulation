package ch.ethz.khammash.hybridstochasticsimulation.simulators;

import java.util.LinkedList;
import java.util.List;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.util.FastMath;

import ch.ethz.khammash.hybridstochasticsimulation.models.AdaptiveMSHRNModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;

abstract public class StepPDMPSimulator extends AbstractSimulator<PDMPModel> {

	protected RandomDataGenerator rdg;
	protected List<TrajectoryRecorder> trajectoryRecorders;
	protected List<TrajectoryRecorder> optionalTrajectoryRecorders;
	protected List<TrajectoryRecorder> simulationInformationTrajectoryRecorders;

	public StepPDMPSimulator(RandomDataGenerator rdg) {
		if (rdg == null)
			rdg = new RandomDataGenerator();
		this.rdg = rdg;
		trajectoryRecorders = new LinkedList<TrajectoryRecorder>();
		optionalTrajectoryRecorders = new LinkedList<TrajectoryRecorder>();
		simulationInformationTrajectoryRecorders = new LinkedList<TrajectoryRecorder>();
	}

	protected double computeUnitJumpTime() {
		return - FastMath.log(rdg.nextUniform(0.0,  1.0));
	}

	@Override
	public void addTrajectoryRecorder(TrajectoryRecorder tr) {
		trajectoryRecorders.add(tr);
	}

	@Override
	public void removeTrajectoryRecorder(TrajectoryRecorder tr) {
		trajectoryRecorders.remove(tr);
	}

	@Override
	public void clearTrajectoryRecorders() {
		trajectoryRecorders.clear();
	}

	public void addOptionalTrajectoryRecorder(TrajectoryRecorder tr) {
		optionalTrajectoryRecorders.add(tr);
	}

	public void removeOptionalTrajectoryRecorder(TrajectoryRecorder tr) {
		optionalTrajectoryRecorders.remove(tr);
	}

	public void clearOptionalTrajectoryRecorders() {
		optionalTrajectoryRecorders.clear();
	}

	public void addSimulationInformationTrajectoryRecorder(TrajectoryRecorder tr) {
		simulationInformationTrajectoryRecorders.add(tr);
	}

	public void removeSimulationInformationTrajectoryRecorder(TrajectoryRecorder tr) {
		simulationInformationTrajectoryRecorders.remove(tr);
	}

	public void clearSimulationInformationTrajectoryRecorders() {
		simulationInformationTrajectoryRecorders.clear();
	}

	protected void beginRecording(AdaptiveMSHRNModel model, PDMPSimulationInformation simInfo, double t0, double[] x0, double t1) {
		double[] primaryState = model.computePrimaryState(t0, x0);
    	for (TrajectoryRecorder tr : trajectoryRecorders)
    		tr.beginRecording(t0, primaryState, t1);
		if (model.hasOptionalState()) {
			double[] optionalState = model.computeOptionalState(t0, x0);
			for (TrajectoryRecorder tr : optionalTrajectoryRecorders)
				tr.beginRecording(t0, optionalState, t1);
		}
		if (simInfo != null) {
			double[] informationState = simInfo.computeInformationState();
			for (TrajectoryRecorder tr : simulationInformationTrajectoryRecorders)
				tr.beginRecording(t0, informationState, t1);
		}
	}

	protected void record(AdaptiveMSHRNModel model, PDMPSimulationInformation simInfo, double t, double[] x) {
		double[] primaryState = model.computePrimaryState(t, x);
		for (TrajectoryRecorder tr : trajectoryRecorders)
			tr.record(t, primaryState);
		if (model.hasOptionalState()) {
			double[] optionalState = model.computeOptionalState(t, x);
			for (TrajectoryRecorder tr : optionalTrajectoryRecorders)
				tr.record(t, optionalState);
		}
		if (simInfo != null) {
			double[] informationState = simInfo.computeInformationState();
			for (TrajectoryRecorder tr : simulationInformationTrajectoryRecorders)
				tr.record(t, informationState);
		}
	}

	protected void endRecording(AdaptiveMSHRNModel model, PDMPSimulationInformation simInfo, double t1, double[] x1) {
		double[] primaryState = model.computePrimaryState(t1, x1);
    	for (TrajectoryRecorder tr : trajectoryRecorders)
    		tr.endRecording(primaryState);
		if (model.hasOptionalState()) {
			double[] optionalState = model.computeOptionalState(t1, x1);
			for (TrajectoryRecorder tr : optionalTrajectoryRecorders)
				tr.endRecording(optionalState);
		}
		if (simInfo != null) {
			double[] informationState = simInfo.computeInformationState();
			for (TrajectoryRecorder tr : simulationInformationTrajectoryRecorders)
				tr.endRecording(informationState);
		}
	}

}
