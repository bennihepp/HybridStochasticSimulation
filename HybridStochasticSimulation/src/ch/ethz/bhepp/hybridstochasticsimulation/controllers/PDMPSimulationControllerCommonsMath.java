package ch.ethz.bhepp.hybridstochasticsimulation.controllers;
//package ch.ethz.khammash.hybridstochasticsimulation.controllers;
//
//import java.util.concurrent.ExecutorService;
//
//import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
//import ch.ethz.khammash.hybridstochasticsimulation.providers.DormandPrince853IntegratorProvider;
//import ch.ethz.khammash.hybridstochasticsimulation.providers.IntegratorProvider;
//import ch.ethz.khammash.hybridstochasticsimulation.providers.PDMPSimulatorProviderCommonsMath;
//import ch.ethz.khammash.hybridstochasticsimulation.providers.UnivariateSolverProvider;
//
//
//public class PDMPSimulationControllerCommonsMath
//		extends AbstractSimulationController<PDMPModel> {
//
//    public PDMPSimulationControllerCommonsMath() {
//    	super();
//		construct();
//    }
//
//	public PDMPSimulationControllerCommonsMath(int numOfThreads) {
//		super(numOfThreads);
//		construct();
//    }
//
//	public PDMPSimulationControllerCommonsMath(ExecutorService executor) {
//		super(executor);
//		construct();
//    }
//
//	final private void construct() {
//		usePDMPSimulatorFactory();
//	}
//
//	final public void usePDMPSimulatorFactory() {
//		usePDMPSimulatorFactory(new DormandPrince853IntegratorProvider());
//	}
//
//	final public void usePDMPSimulatorFactory(IntegratorProvider integratorFactory) {
//		setSimulatorProvider(new PDMPSimulatorProviderCommonsMath(integratorFactory));
//	}
//
//	final public void usePDMPSimulatorFactory(UnivariateSolverProvider univariateSolverFactory) {
//		usePDMPSimulatorFactory(new DormandPrince853IntegratorProvider(), univariateSolverFactory);
//	}
//
//	final public void usePDMPSimulatorFactory(IntegratorProvider integratorFactory, UnivariateSolverProvider univariateSolverFactory) {
//		setSimulatorProvider(new PDMPSimulatorProviderCommonsMath(integratorFactory, univariateSolverFactory));
//	}
//
//}
