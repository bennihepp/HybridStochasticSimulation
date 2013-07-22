package ch.ethz.khammash.hybridstochasticsimulation.batch.guiceproviders;

import org.apache.commons.configuration.HierarchicalConfiguration;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.controllers.PDMPSimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.controllers.SimulationController;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.providers.ObjProvider;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;

import com.google.inject.Inject;

public class PDMPSimulationControllerProvider extends AbstractSimulationControllerProvider<PDMPModel> {

	private ObjProvider<Simulator<PDMPModel>> simulatorProvider;

	@Inject
	public PDMPSimulationControllerProvider(HierarchicalConfiguration config, ObjProvider<Simulator<PDMPModel>> simulatorProvider,
			ObjProvider<RandomDataGenerator> rdgProvider) {
		super(config, rdgProvider);
		this.simulatorProvider = simulatorProvider;
	}

	@Override
	protected SimulationController<PDMPModel> getController(int numOfThreads) {
		PDMPSimulationController simCtrl = new PDMPSimulationController(numOfThreads);
		simCtrl.setSimulatorProvider(simulatorProvider);
		return simCtrl;
	}

}