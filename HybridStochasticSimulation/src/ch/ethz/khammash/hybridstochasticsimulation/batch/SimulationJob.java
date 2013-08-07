package ch.ethz.khammash.hybridstochasticsimulation.batch;

import java.io.IOException;

import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;


public interface SimulationJob {

	public static class OutputAlreadyWrittenException extends RuntimeException {

		private static final long serialVersionUID = -1303336105879923868L;

		public OutputAlreadyWrittenException(String msg) {
			super(msg);
		}

	}

	public enum Type {
		TRAJECTORY, DISTRIBUTION
	}

//	List<SimulationOutput> getOutputs();
//
//	SimulationController<T> getSimulationController();
//
//	T createModel();
//
//	ModelFactory<T> getModelFactory();
//
//	FiniteTrajectoryRecorder createTrajectory();
//
//	FiniteTrajectoryRecorderFactory getTrajectoryFactory();
//
//	double gett0();
//
//	double gett1();
//
//	double[] getx0();
//
	int getRuns();
//
//	double[] getPlotScales();
//
//	String[] getLabels();
//
//	Type getSimulationType();

	void runJob();

	FiniteTrajectory runSingleSimulation();

	void writeOutputs() throws IOException;

	void addSimulationResult(FiniteTrajectory tr);

}
