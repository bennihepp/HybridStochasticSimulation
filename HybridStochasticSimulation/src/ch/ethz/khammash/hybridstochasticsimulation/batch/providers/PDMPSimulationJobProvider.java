package ch.ethz.khammash.hybridstochasticsimulation.batch.providers;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.batch.DefaultSimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob;
import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob.Type;
import ch.ethz.khammash.hybridstochasticsimulation.controllers.PDMPSimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.factories.FiniteTrajectoryRecorderFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.ModelFactory;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;

import com.google.inject.Inject;

public class PDMPSimulationJobProvider extends AbstractSimulationJobProvider {

	private ModelFactory<PDMPModel> modelFactory;
	private FiniteTrajectoryRecorderFactory trajectoryRecorderFactory;
	private PDMPSimulationController simulationController;

	@Inject
	public PDMPSimulationJobProvider(HierarchicalConfiguration config,
			ModelFactory<PDMPModel> modelFactory,
			FiniteTrajectoryRecorderFactory trajectoryRecorderFactory,
			PDMPSimulationController simulationController) {
		super(config);
		this.modelFactory = modelFactory;
		this.trajectoryRecorderFactory = trajectoryRecorderFactory;
		this.simulationController = simulationController;
	}

	@Override
	protected SimulationJob get(double t0, double t1, double[] x0, int runs, Type simulationType) {
		return new DefaultSimulationJob<>(modelFactory, trajectoryRecorderFactory, simulationController, t0, t1, x0, runs, simulationType);
	}

}
