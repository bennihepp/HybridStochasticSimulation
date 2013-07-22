package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import java.util.concurrent.ExecutorService;

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

	// TODO
	final private void construct() {
//		usePDMPSimulatorFactory();
	}

//	final public void usePDMPSimulatorFactory() {
//		usePDMPSimulatorFactory(new CVodeSolverProvider(1e-3, 1e-3));
//	}
//
//	final public void usePDMPSimulatorFactory(ObjProvider<Solver> solverFactory) {
//		setSimulatorProvider(new PDMPSimulatorProvider(solverFactory));
//	}

}
