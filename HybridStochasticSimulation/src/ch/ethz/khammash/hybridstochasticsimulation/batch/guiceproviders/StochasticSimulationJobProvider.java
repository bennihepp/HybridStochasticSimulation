package ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders;

import javax.inject.Provider;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.batch.DefaultSimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob.Type;
import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationOutput;
import ch.ethz.khammash.hybridstochasticsimulation.controllers.SimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.providers.ObjProvider;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;

import com.google.inject.Inject;

public class StochasticSimulationJobProvider extends AbstractSimulationJobProvider {

	private ObjProvider<StochasticReactionNetworkModel> modelProvider;
	private ObjProvider<FiniteTrajectoryRecorder> trajectoryRecorderProvider;
	private Provider<SimulationController<StochasticReactionNetworkModel>> simulationControllerProvider;

	@Inject
	public StochasticSimulationJobProvider(HierarchicalConfiguration config,
			ObjProvider<StochasticReactionNetworkModel> modelProvider,
			ObjProvider<FiniteTrajectoryRecorder> trajectoryRecorderProvider,
			Provider<SimulationController<StochasticReactionNetworkModel>> simulationControllerProvider,
			Provider<SimulationOutput> outputProvider) {
		super(config, outputProvider);
		this.modelProvider = modelProvider;
		this.trajectoryRecorderProvider = trajectoryRecorderProvider;
		this.simulationControllerProvider = simulationControllerProvider;
	}

	@Override
	protected DefaultSimulationJob<?> getSimulationJob(double t0, double t1, double[] x0, int runs, Type simulationType) {
		SimulationController<StochasticReactionNetworkModel> simulationController = simulationControllerProvider.get();
		return new DefaultSimulationJob<>(modelProvider, trajectoryRecorderProvider, simulationController, t0, t1, x0, runs, simulationType);
	}

}

