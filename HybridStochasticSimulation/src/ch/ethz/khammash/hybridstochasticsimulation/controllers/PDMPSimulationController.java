package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import java.util.concurrent.ExecutorService;

import ch.ethz.khammash.hybridstochasticsimulation.factories.CVodeSolverFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.PDMPSimulatorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.SolverFactory;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;


public class PDMPSimulationController extends AbstractSimulationController<PDMPModel> {

    public PDMPSimulationController() {
    	super();
		construct();
    }

	public PDMPSimulationController(int numOfThreads) {
		super(numOfThreads);
		construct();
    }

	public PDMPSimulationController(ExecutorService executor) {
		super(executor);
		construct();
    }

	final private void construct() {
		usePDMPSimulatorFactory();
	}

	final public void usePDMPSimulatorFactory() {
		usePDMPSimulatorFactory(new CVodeSolverFactory(1e-3, 1e-3));
	}

	final public void usePDMPSimulatorFactory(SolverFactory solverFactory) {
		setSimulatorFactory(new PDMPSimulatorFactory(solverFactory));
	}

}
