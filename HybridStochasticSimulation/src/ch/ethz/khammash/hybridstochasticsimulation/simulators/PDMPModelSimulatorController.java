package ch.ethz.khammash.hybridstochasticsimulation.simulators;

import java.util.concurrent.Callable;
import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.ode.AbstractIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.StatisticalSummary;
import org.apache.commons.math3.stat.descriptive.SynchronizedSummaryStatistics;

import ch.ethz.khammash.hybridstochasticsimulation.models.FinitePDMPModelTrajectory;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModelTrajectory;


public class PDMPModelSimulatorController {

	public static interface PDMPModelFactory {
		public PDMPModel createModel();
	}

	public static interface FinitePDMPModelTrajectoryFactory {
		public FinitePDMPModelTrajectory createModelTrajectory();
	}

	public static interface PDMPModelAndFiniteTrajectoryFactory {
		public void createNextModelAndTrajectory();
		public PDMPModel getModel();
		public FinitePDMPModelTrajectory getModelTrajectory();
	}

	public static interface IntegratorFactory {
		public AbstractIntegrator createIntegrator();
	}

	public static class DefaultIntegratorFactory implements IntegratorFactory {

		private double minStep;
		private double maxStep;
        private double scalAbsoluteTolerance;
        private double scalRelativeTolerance;
        private int maxEvaluations;

		public DefaultIntegratorFactory() {
			this(1.0e-10, 100.0, 1.0e-10, 1.0e-10, Integer.MAX_VALUE);
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
		public AbstractIntegrator createIntegrator() {
			DormandPrince853Integrator integrator = new DormandPrince853Integrator(getMinStep(), getMaxStep(), getScalAbsoluteTolerance(), getScalRelativeTolerance());
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

	private class SimulationWorker implements Callable<PDMPModelTrajectory> {

		private PDMPModel model;
		private double t0;
		private double[] x0;
		private double t1;
		private PDMPModelSimulator simulator;

		public SimulationWorker(AbstractIntegrator integrator,
				RandomDataGenerator rng, PDMPModel model, double t0,
				double[] x0, double t1) {
			this.model = model;
			this.t0 = t0;
			this.x0 = x0;
			this.t1 = t1;
			simulator = new PDMPModelSimulator(integrator, rng);
		}

		public PDMPModelTrajectory simulate() {
			double[] x1 = new double[x0.length];
			PDMPModelTrajectory cm = new PDMPModelTrajectory();
			simulator.addStepHandler(cm);
			simulator.addReactionHandler(cm);
			simulator.simulate(model, t0, x0, t1, x1);
			simulator.clearReactionHandlers();
			simulator.clearStepHandlers();
			return cm;
		}

		@Override
		public PDMPModelTrajectory call() {
			return simulate();
		}
	}

	private class FiniteSimulationWorker implements Callable<FinitePDMPModelTrajectory> {

		private FinitePDMPModelTrajectory mt;
		private PDMPModel model;
		private double t0;
		private double[] x0;
		private double t1;
		private PDMPModelSimulator simulator;

		public FiniteSimulationWorker(FinitePDMPModelTrajectory mt,
				AbstractIntegrator integrator, RandomDataGenerator rdg,
				PDMPModel model, double t0, double[] x0, double t1) {
			this.mt = mt;
			this.model = model;
			this.t0 = t0;
			this.x0 = x0;
			this.t1 = t1;
			simulator = new PDMPModelSimulator(integrator, rdg);
		}

		public FinitePDMPModelTrajectory simulate() {
			double[] x1 = new double[x0.length];
			simulator.addStepHandler(mt);
			simulator.addReactionHandler(mt);
			simulator.simulate(model, t0, x0, t1, x1);
			simulator.clearReactionHandlers();
			simulator.clearStepHandlers();
			return mt;
		}

		@Override
		public FinitePDMPModelTrajectory call() {
			return simulate();
		}
	}

	private class TrajectoryDistributionSimulator extends SimulationWorker implements Runnable {
		private double[] tSeries;
		private SynchronizedSummaryStatistics[][] xStatistics;

		public TrajectoryDistributionSimulator(AbstractIntegrator integrator,
				RandomDataGenerator rng, PDMPModel model,
				SynchronizedSummaryStatistics[][] xStatistics,
				double[] tSeries, double[] x0) {
			super(integrator, rng, model, tSeries[0], x0, tSeries[tSeries.length - 1]);
			this.xStatistics = xStatistics;
			this.tSeries = tSeries;
		}

		@Override
		public void run() {
			PDMPModelTrajectory cm = simulate();
			for (int i = 0; i < tSeries.length; i++) {
				cm.setInterpolatedTime(tSeries[i]);
				double[] x = cm.getInterpolatedState();
				for (int s = 0; s < x.length; s++) {
					xStatistics[i][s].addValue(x[s]);
				}
			}
		}
	}

	private class FiniteTrajectoryDistributionSimulator extends FiniteSimulationWorker implements Runnable {
		private SynchronizedSummaryStatistics[][] xStatistics;

		public FiniteTrajectoryDistributionSimulator(FinitePDMPModelTrajectory mt,
				AbstractIntegrator integrator, RandomDataGenerator rdg,
				PDMPModel model, SynchronizedSummaryStatistics[][] xStatistics,
				double t0, double[] x0, double t1) {
			super(mt, integrator, rdg, model, t0, x0, t1);
			this.xStatistics = xStatistics;
		}

		@Override
		public void run() {
			FinitePDMPModelTrajectory mt = simulate();
			double[][] xSeries = mt.getxSeries();
			for (int s = 0; s < xStatistics.length; s++)
				for (int i = 0; i < xStatistics[s].length; i++)
					xStatistics[s][i].addValue(xSeries[s][i]);
		}
	}

	public final static int DEFAULT_NUMBER_OF_THREADS = 4;

	private IntegratorFactory integratorFactory;
	private ExecutorService executor;
	private RandomGenerator rng;

	public PDMPModelSimulatorController() {
		this(DEFAULT_NUMBER_OF_THREADS);
	}

	public PDMPModelSimulatorController(int numOfThreads) {
		this(Executors.newFixedThreadPool(numOfThreads));
	}

	public PDMPModelSimulatorController(ExecutorService executor) {
		integratorFactory = new DefaultIntegratorFactory();
		this.executor = executor;
	}

	public void setIntegratorFactory(IntegratorFactory integratorFactory) {
		this.integratorFactory = integratorFactory;
	}

	public void setRandomGenerator(RandomGenerator rng) {
		this.rng = rng;
	}

	private AbstractIntegrator createIntegrator() {
		return integratorFactory.createIntegrator();
	}

// TODO
//	private AbstractIntegrator[] createIntegrators(int numOfGenerators) {
//		AbstractIntegrator[] result = new AbstractIntegrator[numOfGenerators];
//		for (int i = 0; i < result.length; i++)
//			result[i] = createIntegrator();
//		return result;
//	}

	private RandomDataGenerator createRandomDataGenerator() {
		if (rng == null)
			rng = new MersenneTwister();
		RandomDataGenerator rdg = new RandomDataGenerator();
		rdg.reSeed(rng.nextLong());
		return rdg;
	}

// TODO
//	private RandomDataGenerator[] createRandomDataGenerators(int numOfGenerators) {
//		RandomDataGenerator[] result = new RandomDataGenerator[numOfGenerators];
//		for (int i = 0; i < result.length; i++)
//			result[i] = createRandomDataGenerator();
//		return result;
//	}

	public PDMPModelTrajectory simulateTrajectory(PDMPModel model, double t0, double[] x0, double t1) {
		RandomDataGenerator rdg = createRandomDataGenerator();
		AbstractIntegrator integrator = createIntegrator();
		SimulationWorker sw = new SimulationWorker(integrator, rdg, model, t0, x0, t1);
		return sw.simulate();
	}

	public double[][] simulateTrajectory(PDMPModel model, double[] tSeries, double[] x0) {
		FinitePDMPModelTrajectory mt = new FinitePDMPModelTrajectory(tSeries);
		simulateTrajectory(model, mt, tSeries[0], x0, tSeries[tSeries.length - 1]);
		return mt.getxSeries();
	}

	public void simulateTrajectory(PDMPModel model, FinitePDMPModelTrajectory mt, double t0, double[] x0, double t1) {
		RandomDataGenerator rdg = createRandomDataGenerator();
		AbstractIntegrator integrator = createIntegrator();
		FiniteSimulationWorker sw = new FiniteSimulationWorker(mt, integrator, rdg, model, t0, x0, t1);
		sw.simulate();
	}

	public StatisticalSummary[][] simulateTrajectoryDistribution(int runs,
			PDMPModelFactory modelFactory, double[] tSeries, double[] x0)
			throws InterruptedException, CancellationException, ExecutionException {
		SynchronizedSummaryStatistics[][] xSeriesStatistics = new SynchronizedSummaryStatistics[tSeries.length][x0.length];
		for (int i = 0; i < tSeries.length; i++)
			for (int s = 0; s < x0.length; s++)
				xSeriesStatistics[i][s] = new SynchronizedSummaryStatistics();
		for (int k = 0; k < runs; k++) {
			AbstractIntegrator integrator = createIntegrator();
			RandomDataGenerator rdg = createRandomDataGenerator();
			// PDMPModel has an internal state so we create a copy for each run
			PDMPModel model = modelFactory.createModel();
			Runnable r = new TrajectoryDistributionSimulator(integrator, rdg, model, xSeriesStatistics, tSeries, x0);
			executor.execute(r);
		}
		executor.shutdown();
		try {
			executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
		} catch (InterruptedException e) {
			executor.shutdownNow();
			Thread.currentThread().interrupt();
		}
		return xSeriesStatistics;
	}

	public StatisticalSummary[][] simulateFiniteTrajectoryDistribution(int runs,
			PDMPModelFactory modelFactory,
			FinitePDMPModelTrajectoryFactory mtFactory, double t0, double[] x0,
			double t1) throws InterruptedException, CancellationException,
			ExecutionException {
		SynchronizedSummaryStatistics[][] xSeriesStatistics = null;
		// We also want to have different random number sequences for each run
		for (int k = 0; k < runs; k++) {
			AbstractIntegrator integrator = createIntegrator();
			RandomDataGenerator rdg = createRandomDataGenerator();
			// PDMPModel has an internal state so we create a copy for each run
			PDMPModel model = modelFactory.createModel();
			FinitePDMPModelTrajectory mt = mtFactory.createModelTrajectory();
			if (xSeriesStatistics == null) {
				xSeriesStatistics = new SynchronizedSummaryStatistics[x0.length][mt.gettSeries().length];
				for (int s = 0; s < xSeriesStatistics.length; s++)
					for (int i = 0; i < xSeriesStatistics[s].length; i++)
						xSeriesStatistics[s][i] = new SynchronizedSummaryStatistics();
			}
			Runnable r = new FiniteTrajectoryDistributionSimulator(mt,
					integrator, rdg, model, xSeriesStatistics, t0, x0.clone(), t1);
			executor.execute(r);
		}
		executor.shutdown();
		try {
			executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
		} catch (InterruptedException e) {
			executor.shutdownNow();
			Thread.currentThread().interrupt();
		}
		return xSeriesStatistics;
	}

}
