package ch.ethz.khammash.hybridstochasticsimulation.batch.providers;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.factories.PDMPSimulatorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.SimulatorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.SolverFactory;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;

import com.google.inject.Inject;

public class PDMPSimulatorFactoryProvider extends AbstractProvider<SimulatorFactory<PDMPModel>> {

	private SolverFactory solverFactory;

	@Inject
	public PDMPSimulatorFactoryProvider(HierarchicalConfiguration config, SolverFactory solverFactory) {
		super(config, "SimulationParameters");
		this.solverFactory = solverFactory;
	}

	@Override
	public SimulatorFactory<PDMPModel> get() {
		return new PDMPSimulatorFactory(solverFactory);
	}

}
