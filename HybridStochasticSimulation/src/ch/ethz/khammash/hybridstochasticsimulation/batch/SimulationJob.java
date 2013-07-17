package ch.ethz.khammash.hybridstochasticsimulation.batch;


public interface SimulationJob {

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
//	int getRuns();
//
//	double[] getPlotScales();
//
//	String[] getLabels();
//
//	Type getSimulationType();

	void runJob();

}
