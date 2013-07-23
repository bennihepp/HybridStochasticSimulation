package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import java.util.concurrent.ExecutorService;

import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.providers.RandomDataGeneratorProvider;
import ch.ethz.khammash.hybridstochasticsimulation.providers.StochasticSimulatorProvider;


public class StochasticSimulationController extends AbstractSimulationController<StochasticReactionNetworkModel> {

    public StochasticSimulationController() {
    	super();
		construct();
    }

	public StochasticSimulationController(int numOfThreads) {
		super(numOfThreads);
		construct();
    }

	public StochasticSimulationController(ExecutorService executor) {
		super(executor);
		construct();
    }

	final private void construct() {
//		setSimulatorProvider(new StochasticSimulatorProvider());
	}

	public void setDefaultSimulatorProvider(RandomDataGeneratorProvider rdgProvider) {
		setRandomDataGeneratorProvider(rdgProvider);
		setSimulatorProvider(new StochasticSimulatorProvider(rdgProvider));
		
	}

}
