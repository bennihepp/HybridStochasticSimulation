package ch.ethz.khammash.hybridstochasticsimulation.simulators;

import java.util.concurrent.Callable;

import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;

import ch.ethz.khammash.hybridstochasticsimulation.models.DeterministicModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.FiniteDeterministicModelTrajectory;
import ch.ethz.khammash.hybridstochasticsimulation.models.ModelTrajectory;


public class DeterministicModelSimulatorController {

	public static interface IntegratorFactory {
		public FirstOrderIntegrator createIntegrator();
	}

	public static class DefaultIntegratorFactory implements IntegratorFactory {

		private double minStep;
		private double maxStep;
        private double scalAbsoluteTolerance;
        private double scalRelativeTolerance;
        private int maxEvaluations;

		public DefaultIntegratorFactory() {
			this(1.0e-10, 100.0, 1.0e-6, 1.0e-3, Integer.MAX_VALUE);
		}

		public DefaultIntegratorFactory(double minStep, double maxStep,
                double scalAbsoluteTolerance, double scalRelativeTolerance,
                int maxEvaluations) {
			this.setMinStep(minStep);
			this.setMaxStep(maxStep);
			this.setScalAbsoluteTolerance(scalAbsoluteTolerance);
			this.setScalRelativeTolerance(scalRelativeTolerance);
			this.setMaxEvaluations(maxEvaluations);
		}

		@Override
		public FirstOrderIntegrator createIntegrator() {
			DormandPrince853Integrator integrator = new DormandPrince853Integrator(
					getMinStep(), getMaxStep(), getScalAbsoluteTolerance(),
					getScalRelativeTolerance());
			integrator.setMaxEvaluations(maxEvaluations);
			return integrator;
		}

		public double getMinStep() {
			return minStep;
		}

		public void setMinStep(double minStep) {
			this.minStep = minStep;
		}

		public double getMaxStep() {
			return maxStep;
		}

		public void setMaxStep(double maxStep) {
			this.maxStep = maxStep;
		}

		public double getScalAbsoluteTolerance() {
			return scalAbsoluteTolerance;
		}

		public void setScalAbsoluteTolerance(double scalAbsoluteTolerance) {
			this.scalAbsoluteTolerance = scalAbsoluteTolerance;
		}

		public double getScalRelativeTolerance() {
			return scalRelativeTolerance;
		}

		public void setScalRelativeTolerance(double scalRelativeTolerance) {
			this.scalRelativeTolerance = scalRelativeTolerance;
		}

		public int getMaxEvaluations() {
			return maxEvaluations;
		}

		public void setMaxEvaluations(int maxEvaluations) {
			this.maxEvaluations = maxEvaluations;
		}
	}

	private class SimulationWorker implements Callable<ModelTrajectory> {

		private FiniteDeterministicModelTrajectory mt;
		private DeterministicModel model;
		private double t0;
		private double[] x0;
		private double t1;
		private DeterministicModelSimulator simulator;

		public SimulationWorker(FiniteDeterministicModelTrajectory mt,
				FirstOrderIntegrator integrator, DeterministicModel model,
				double t0, double[] x0, double t1) {
			this.mt = mt;
			this.model = model;
			this.t0 = t0;
			this.x0 = x0;
			this.t1 = t1;
			simulator = new DeterministicModelSimulator(integrator);
		}

		public FiniteDeterministicModelTrajectory simulate() {
			double[] x1 = new double[x0.length];
			simulator.addStepHandler(mt);
			simulator.simulate(model, t0, x0, t1, x1);
			simulator.clearStepHandlers();
			return mt;
		}

		@Override
		public FiniteDeterministicModelTrajectory call() {
			return simulate();
		}
	}

	public final static int DEFAULT_NUMBER_OF_THREADS = 4;

	private IntegratorFactory integratorFactory;

	public DeterministicModelSimulatorController() {
		integratorFactory = new DefaultIntegratorFactory();
	}

	public void setIntegratorFactory(IntegratorFactory integratorFactory) {
		this.integratorFactory = integratorFactory;
	}

	private FirstOrderIntegrator createIntegrator() {
		return integratorFactory.createIntegrator();
	}

	public FiniteDeterministicModelTrajectory simulateTrajectory(DeterministicModel model, double[] tSeries, double[] x0) {
		FiniteDeterministicModelTrajectory mt = new FiniteDeterministicModelTrajectory(tSeries);
		FirstOrderIntegrator integrator = createIntegrator();
		SimulationWorker sw = new SimulationWorker(mt, integrator, model, tSeries[0], x0, tSeries[tSeries.length - 1]);
		return sw.simulate();
	}

}
