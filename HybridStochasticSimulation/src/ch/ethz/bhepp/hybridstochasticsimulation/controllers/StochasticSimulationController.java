package ch.ethz.bhepp.hybridstochasticsimulation.controllers;

import java.util.concurrent.ExecutorService;

import ch.ethz.bhepp.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.bhepp.hybridstochasticsimulation.providers.ObjProvider;
import ch.ethz.bhepp.hybridstochasticsimulation.simulators.Simulator;


public class StochasticSimulationController extends AbstractSimulationController<StochasticReactionNetworkModel> {

    public StochasticSimulationController(ObjProvider<? extends Simulator<StochasticReactionNetworkModel>> simulatorProvider) {
    	super(simulatorProvider);
    }

	public StochasticSimulationController(ObjProvider<? extends Simulator<StochasticReactionNetworkModel>> simulatorProvider,
			int numOfThreads) {
		super(simulatorProvider, numOfThreads);
    }

	public StochasticSimulationController(ObjProvider<? extends Simulator<StochasticReactionNetworkModel>> simulatorProvider,
			ExecutorService executor) {
		super(simulatorProvider, executor);
    }

}
