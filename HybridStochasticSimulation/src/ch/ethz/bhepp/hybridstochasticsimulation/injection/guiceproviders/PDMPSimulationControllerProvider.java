package ch.ethz.bhepp.hybridstochasticsimulation.injection.guiceproviders;

import org.apache.commons.configuration.HierarchicalConfiguration;

import ch.ethz.bhepp.hybridstochasticsimulation.controllers.PDMPSimulationController;
import ch.ethz.bhepp.hybridstochasticsimulation.controllers.SimulationController;
import ch.ethz.bhepp.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.bhepp.hybridstochasticsimulation.providers.ObjProvider;
import ch.ethz.bhepp.hybridstochasticsimulation.simulators.Simulator;

import com.google.inject.Inject;

public class PDMPSimulationControllerProvider extends AbstractSimulationControllerProvider<PDMPModel> {

	private ObjProvider<Simulator<PDMPModel>> simulatorProvider;

	@Inject
	public PDMPSimulationControllerProvider(HierarchicalConfiguration config, ObjProvider<Simulator<PDMPModel>> simulatorProvider) {
		super(config);
		this.simulatorProvider = simulatorProvider;
	}

	@Override
	protected SimulationController<PDMPModel> getController(int numOfThreads) {
		PDMPSimulationController simCtrl = new PDMPSimulationController(simulatorProvider, numOfThreads);
		return simCtrl;
	}

}
