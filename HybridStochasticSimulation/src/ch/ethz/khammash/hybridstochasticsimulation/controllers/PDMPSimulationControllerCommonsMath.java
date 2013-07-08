package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import java.util.concurrent.ExecutorService;

import ch.ethz.khammash.hybridstochasticsimulation.factories.DormandPrince853IntegratorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.IntegratorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.PDMPSimulatorFactoryCommonsMath;
import ch.ethz.khammash.hybridstochasticsimulation.factories.UnivariateSolverFactory;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;


public class PDMPSimulationControllerCommonsMath
		extends AbstractSimulationController<PDMPModel> {

    public PDMPSimulationControllerCommonsMath() {
    	super();
		construct();
    }

	public PDMPSimulationControllerCommonsMath(int numOfThreads) {
		super(numOfThreads);
		construct();
    }

	public PDMPSimulationControllerCommonsMath(ExecutorService executor) {
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
		setSimulatorFactory(new PDMPSimulatorFactoryCommonsMath(integratorFactory));
	}

	final public void usePDMPSimulatorFactory(UnivariateSolverFactory univariateSolverFactory) {
		usePDMPSimulatorFactory(new DormandPrince853IntegratorFactory(), univariateSolverFactory);
	}

	final public void usePDMPSimulatorFactory(IntegratorFactory integratorFactory, UnivariateSolverFactory univariateSolverFactory) {
		setSimulatorFactory(new PDMPSimulatorFactoryCommonsMath(integratorFactory, univariateSolverFactory));
	}

}
