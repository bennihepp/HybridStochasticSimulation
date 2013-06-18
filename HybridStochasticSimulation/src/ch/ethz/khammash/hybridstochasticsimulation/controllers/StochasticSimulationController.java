package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import java.util.concurrent.ExecutorService;

import ch.ethz.khammash.hybridstochasticsimulation.factories.StochasticSimulatorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;


public class StochasticSimulationController<T extends StochasticReactionNetworkModel>
		extends AbstractSimulationController<T, TrajectoryRecorder<T>> {

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
		setSimulatorFactory(new StochasticSimulatorFactory<T>());
	}

}
