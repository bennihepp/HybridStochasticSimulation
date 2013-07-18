package ch.ethz.khammash.hybridstochasticsimulation.batch.providers;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.khammash.hybridstochasticsimulation.controllers.PDMPSimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.controllers.SimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.factories.RandomDataGeneratorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.SimulatorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;

import com.google.inject.Inject;

public class PDMPSimulationControllerProvider extends AbstractSimulationControllerProvider<PDMPModel> {

	private SimulatorFactory<PDMPModel> simulatorFactory;

	@Inject
	public PDMPSimulationControllerProvider(HierarchicalConfiguration config, SimulatorFactory<PDMPModel> simulatorFactory,
			RandomDataGeneratorFactory rdgFactory) {
		super(config, rdgFactory);
		this.simulatorFactory = simulatorFactory;
	}

	@Override
	protected SimulationController<PDMPModel> getController(int numOfThreads) {
		PDMPSimulationController simCtrl = new PDMPSimulationController(numOfThreads);
		simCtrl.setSimulatorFactory(simulatorFactory);
		return simCtrl;
	}

}
