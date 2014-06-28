package ch.ethz.bhepp.hybridstochasticsimulation.controllers;

import java.util.concurrent.ExecutorService;

import ch.ethz.bhepp.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.bhepp.hybridstochasticsimulation.providers.ObjProvider;
import ch.ethz.bhepp.hybridstochasticsimulation.simulators.Simulator;


public class PDMPSimulationController extends AbstractSimulationController<PDMPModel> {

    public PDMPSimulationController(ObjProvider<? extends Simulator<PDMPModel>> simulatorProvider) {
    	super(simulatorProvider);
    }

	public PDMPSimulationController(ObjProvider<? extends Simulator<PDMPModel>> simulatorProvider, int numOfThreads) {
		super(simulatorProvider, numOfThreads);
    }

	public PDMPSimulationController(ObjProvider<? extends Simulator<PDMPModel>> simulatorProvider, ExecutorService executor) {
		super(simulatorProvider, executor);
    }

}
