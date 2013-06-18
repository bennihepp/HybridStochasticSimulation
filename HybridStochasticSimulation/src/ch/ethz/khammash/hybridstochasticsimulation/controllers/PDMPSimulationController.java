package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import java.util.concurrent.ExecutorService;

import ch.ethz.khammash.hybridstochasticsimulation.factories.DormandPrince853IntegratorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.PDMPSimulatorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.IntegratorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.UnivariateSolverFactory;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.ContinuousTrajectoryRecorder;


public class PDMPSimulationController<T extends PDMPModel>
		extends AbstractSimulationController<T, ContinuousTrajectoryRecorder<T>> {

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
		usePDMPSimulatorFactory(new DormandPrince853IntegratorFactory());
	}

	final public void usePDMPSimulatorFactory(IntegratorFactory integratorFactory) {
		setSimulatorFactory(new PDMPSimulatorFactory<T>(integratorFactory));
	}

	final public void usePDMPSimulatorFactory(UnivariateSolverFactory univariateSolverFactory) {
		usePDMPSimulatorFactory(new DormandPrince853IntegratorFactory(), univariateSolverFactory);
	}

	final public void usePDMPSimulatorFactory(IntegratorFactory integratorFactory, UnivariateSolverFactory univariateSolverFactory) {
		setSimulatorFactory(new PDMPSimulatorFactory<T>(integratorFactory, univariateSolverFactory));
	}

}
