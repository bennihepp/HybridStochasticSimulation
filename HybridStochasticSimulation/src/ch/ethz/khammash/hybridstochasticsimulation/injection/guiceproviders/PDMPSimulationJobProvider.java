package ch.ethz.khammash.hybridstochasticsimulation.injection.guiceproviders;



import javax.inject.Provider;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.batch.DefaultSimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob.Type;
import ch.ethz.khammash.hybridstochasticsimulation.controllers.SimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.io.SimulationOutput;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.providers.ObjProvider;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectoryMapper;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;

import com.google.inject.Inject;

public class PDMPSimulationJobProvider extends AbstractSimulationJobProvider {

	private ObjProvider<PDMPModel> modelProvider;
	private ObjProvider<FiniteTrajectoryRecorder> trajectoryRecorderProvider;
	private Provider<SimulationController<PDMPModel>> simulationControllerProvider;

	@Inject
	public PDMPSimulationJobProvider(HierarchicalConfiguration config,
			ObjProvider<PDMPModel> modelProvider,
			ObjProvider<FiniteTrajectoryRecorder> trajectoryRecorderProvider,
			Provider<SimulationController<PDMPModel>> simulationControllerProvider,
			Provider<SimulationOutput> outputProvider,
			Provider<FiniteTrajectoryMapper> mapperProvider) {
		super(config, outputProvider, mapperProvider);
		this.modelProvider = modelProvider;
		this.trajectoryRecorderProvider = trajectoryRecorderProvider;
		this.simulationControllerProvider = simulationControllerProvider;
	}

	@Override
	protected DefaultSimulationJob<?> getSimulationJob(double t0, double t1, double[] x0, int runs, Type simulationType) {
		SimulationController<PDMPModel> simulationController = simulationControllerProvider.get();
		return new DefaultSimulationJob<PDMPModel>(modelProvider, trajectoryRecorderProvider, simulationController, t0, t1, x0, runs, simulationType);
	}

}
