package ch.ethz.khammash.hybridstochasticsimulation.batch.providers;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.batch.DefaultSimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationOutput;
import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob.Type;
import ch.ethz.khammash.hybridstochasticsimulation.controllers.StochasticSimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.factories.FiniteTrajectoryRecorderFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.ModelFactory;
import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;

import com.google.inject.Inject;

public class StochasticSimulationJobProvider extends AbstractSimulationJobProvider {

	private ModelFactory<StochasticReactionNetworkModel> modelFactory;
	private FiniteTrajectoryRecorderFactory trajectoryRecorderFactory;
	private StochasticSimulationController simulationController;

	@Inject
	public StochasticSimulationJobProvider(HierarchicalConfiguration config,
			ModelFactory<StochasticReactionNetworkModel> modelFactory,
			FiniteTrajectoryRecorderFactory trajectoryRecorderFactory,
			StochasticSimulationController simulationController,
			SimulationOutput output) {
		super(config, output);
		this.modelFactory = modelFactory;
		this.trajectoryRecorderFactory = trajectoryRecorderFactory;
		this.simulationController = simulationController;
	}

	@Override
	protected DefaultSimulationJob<?> getSimulationJob(double t0, double t1, double[] x0, int runs, Type simulationType) {
		return new DefaultSimulationJob<>(modelFactory, trajectoryRecorderFactory, simulationController, t0, t1, x0, runs, simulationType);
	}

}

