package ch.ethz.khammash.hybridstochasticsimulation.batch.providers;

import javax.inject.Provider;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.batch.DefaultSimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob.Type;
import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationOutput;
import ch.ethz.khammash.hybridstochasticsimulation.controllers.SimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.factories.FiniteTrajectoryRecorderFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.ModelFactory;
import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;

import com.google.inject.Inject;

public class StochasticSimulationJobProvider extends AbstractSimulationJobProvider {

	private ModelFactory<StochasticReactionNetworkModel> modelFactory;
	private FiniteTrajectoryRecorderFactory trajectoryRecorderFactory;
	private Provider<SimulationController<StochasticReactionNetworkModel>> simulationControllerProvider;

	@Inject
	public StochasticSimulationJobProvider(HierarchicalConfiguration config,
			ModelFactory<StochasticReactionNetworkModel> modelFactory,
			FiniteTrajectoryRecorderFactory trajectoryRecorderFactory,
			Provider<SimulationController<StochasticReactionNetworkModel>> simulationControllerProvider,
			Provider<SimulationOutput> outputProvider) {
		super(config, outputProvider);
		this.modelFactory = modelFactory;
		this.trajectoryRecorderFactory = trajectoryRecorderFactory;
		this.simulationControllerProvider = simulationControllerProvider;
	}

	@Override
	protected DefaultSimulationJob<?> getSimulationJob(double t0, double t1, double[] x0, int runs, Type simulationType) {
		SimulationController<StochasticReactionNetworkModel> simulationController = simulationControllerProvider.get();
		return new DefaultSimulationJob<>(modelFactory, trajectoryRecorderFactory, simulationController, t0, t1, x0, runs, simulationType);
	}

}

