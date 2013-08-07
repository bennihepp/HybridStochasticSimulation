package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import java.util.concurrent.ExecutorService;

import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;


public class PDMPSimulationController extends AbstractSimulationController<PDMPModel> {

    public PDMPSimulationController() {
    	super();
    }

	public PDMPSimulationController(int numOfThreads) {
		super(numOfThreads);
    }

	public PDMPSimulationController(ExecutorService executor) {
		super(executor);
    }

}
