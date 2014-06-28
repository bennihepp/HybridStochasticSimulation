package ch.ethz.bhepp.hybridstochasticsimulation.providers;
//package ch.ethz.khammash.hybridstochasticsimulation.providers;
//
//import org.apache.commons.math3.analysis.solvers.UnivariateSolver;
//import org.apache.commons.math3.ode.AbstractIntegrator;
//import org.apache.commons.math3.random.RandomDataGenerator;
//
//import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
//import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPSimulatorCommonsMath;
//
//public class PDMPSimulatorProviderCommonsMath implements ObjProvider<Simulator<PDMPModel>> {
//
//	private IntegratorProvider integratorFactory;
//	private double ehMaxCheckInterval;
//	private double ehConvergence;
//	private double ehConvergenceFactor;
//	private int ehMaxIterationCount;
//	private UnivariateSolverProvider univariateSolverFactory;
//
//	public PDMPSimulatorProviderCommonsMath(IntegratorProvider integratorFactory) {
//		this(integratorFactory, new BracketingNthOrderBrentSolverProvider());
//	}
//
//	public PDMPSimulatorProviderCommonsMath(IntegratorProvider integratorFactory, UnivariateSolverProvider univariateSolverFactory) {
//		this.integratorFactory = integratorFactory;
//		ehMaxCheckInterval = PDMPSimulatorCommonsMath.DEFAULT_EVENT_HANDLER_MAX_CHECK_INTERVAL;
//		ehConvergence = PDMPSimulatorCommonsMath.DEFAULT_EVENT_HANDLER_CONVERGENCE;
//		ehConvergenceFactor = PDMPSimulatorCommonsMath.DEFAULT_EVENT_HANDLER_CONVERGENCE_FACTOR;
//		ehMaxIterationCount = PDMPSimulatorCommonsMath.DEFAULT_EVENT_HANDLER_MAX_ITERATION_COUNT;
//		this.univariateSolverFactory = univariateSolverFactory;
//	}
//
//	public void setIntegratorFactory(IntegratorProvider integratorFactory) {
//		this.integratorFactory = integratorFactory;
//	}
//
//	public void setUnivariateSolverFactory(UnivariateSolverProvider univariateSolverFactory) {
//		this.univariateSolverFactory = univariateSolverFactory;
//	}
//
//	@Override
//	public PDMPSimulatorCommonsMath createSimulator(RandomDataGenerator rdg) {
//		AbstractIntegrator integrator = integratorFactory.get();
//		UnivariateSolver univariateSolver = univariateSolverFactory.get();
//		PDMPSimulatorCommonsMath sim = new PDMPSimulatorCommonsMath(integrator, univariateSolver, rdg);
//		sim.setEventHandlerConvergence(getEventHandlerConvergence());
//		sim.setEventHandlerConvergenceFactor(getEventHandlerConvergenceFactor());
//		sim.setEventHandlerMaxCheckInterval(getEventHandlerMaxCheckInterval());
//		sim.setEventHandlerMaxIterationCount(getEventHandlerMaxIterationCount());
//		return sim;
//	}
//
//	public double getEventHandlerConvergence() {
//		return ehConvergence;
//	}
//
//	public void setEventHandlerConvergence(double convergence) {
//		this.ehConvergence = convergence;
//		if (!Double.isNaN(convergence))
//			this.ehConvergenceFactor = Double.NaN;
//	}
//
//	public double getEventHandlerConvergenceFactor() {
//		return ehConvergenceFactor;
//	}
//
//	public void setEventHandlerConvergenceFactor(double convergenceFactor) {
//		this.ehConvergenceFactor = convergenceFactor;
//		if (!Double.isNaN(convergenceFactor))
//			this.ehConvergence = Double.NaN;
//	}
//
//	public double getEventHandlerMaxCheckInterval() {
//		return ehMaxCheckInterval;
//	}
//
//	public void setEventHandlerMaxCheckInterval(double maxCheckInterval) {
//		this.ehMaxCheckInterval = maxCheckInterval;
//	}
//
//	public int getEventHandlerMaxIterationCount() {
//		return ehMaxIterationCount;
//	}
//
//	public void setEventHandlerMaxIterationCount(int maxIterationCount) {
//		this.ehMaxIterationCount = maxIterationCount;
//	}
//
//}