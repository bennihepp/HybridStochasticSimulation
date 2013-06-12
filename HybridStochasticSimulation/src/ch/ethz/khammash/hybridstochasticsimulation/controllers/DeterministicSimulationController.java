//package ch.ethz.khammash.hybridstochasticsimulation.controllers;
//
//import java.util.concurrent.Callable;
//
//import org.apache.commons.math3.ode.AbstractIntegrator;
//import org.apache.commons.math3.ode.FirstOrderIntegrator;
//
//import com.google.common.base.Optional;
//
//import ch.ethz.khammash.hybridstochasticsimulation.models.StochasticReactionNetworkModel;
//import ch.ethz.khammash.hybridstochasticsimulation.models.UnaryBinaryDeterministicModel;
//import ch.ethz.khammash.hybridstochasticsimulation.simulators.DeterministicSimulator;
//import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;
//import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteDeterministicTrajectory;
//import ch.ethz.khammash.hybridstochasticsimulation.trajectories.Trajectory;
//import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;
//
//
//public class DeterministicSimulationController extends
//AbstractSimulationController<StochasticReactionNetworkModel, TrajectoryRecorder> {
//
//	private IntegratorFactory integratorFactory;
//
//	public DeterministicSimulationController() {
//		integratorFactory = new DefaultIntegratorFactory();
//	}
//
//	final public void setIntegratorFactory(IntegratorFactory integratorFactory) {
//		this.integratorFactory = integratorFactory;
//	}
//
//	private FirstOrderIntegrator createIntegrator() {
//		return integratorFactory.createIntegrator();
//	}
//
//	public FiniteDeterministicTrajectory simulateTrajectory(UnaryBinaryDeterministicModel model, double[] tSeries, double[] x0) {
//		FiniteDeterministicTrajectory mt = new FiniteDeterministicTrajectory(tSeries);
//		FirstOrderIntegrator integrator = createIntegrator();
//		SimulationWorker sw = new SimulationWorker(mt, integrator, model, tSeries[0], x0, tSeries[tSeries.length - 1]);
//		return sw.simulate();
//	}
//
//}
