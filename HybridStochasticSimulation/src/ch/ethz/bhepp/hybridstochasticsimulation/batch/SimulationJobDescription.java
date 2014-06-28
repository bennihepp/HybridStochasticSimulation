package ch.ethz.bhepp.hybridstochasticsimulation.batch;

import java.util.List;

import ch.ethz.bhepp.hybridstochasticsimulation.controllers.SimulationController;
import ch.ethz.bhepp.hybridstochasticsimulation.io.SimulationOutput;
import ch.ethz.bhepp.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.bhepp.hybridstochasticsimulation.providers.ObjProvider;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;

public interface SimulationJobDescription<T extends ReactionNetworkModel> extends SimulationJob {

	String getName();

	List<SimulationOutput> getOutputs();

	SimulationController<T> getSimulationController();

	T createModel();

	ObjProvider<T> getModelProvider();

	FiniteTrajectoryRecorder createTrajectory();

	ObjProvider<FiniteTrajectoryRecorder> getTrajectoryProvider();

	double gett0();

	double gett1();

	double[] getx0();

//	int getRuns();

	double[] getPlotScales();

	List<String> getLabels();

	Type getSimulationType();

	void runJob() throws InterruptedException;

}
