package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import java.util.concurrent.Callable;
import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.ode.AbstractIntegrator;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.stat.descriptive.StatisticalSummary;
import org.apache.commons.math3.stat.descriptive.SynchronizedSummaryStatistics;

import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPSimulator;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FinitePDMPTrajectory;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.PDMPTrajectory;


public class PDMPSimulatorController {

	private class SimulationWorker implements Callable<PDMPTrajectory> {

		private PDMPModel model;
		private double t0;
		private double[] x0;
		private double t1;
		private PDMPSimulator simulator;

		public SimulationWorker(PDMPSimulator simulator, PDMPModel model, double t0,
				double[] x0, double t1) {
			this.model = model;
			this.t0 = t0;
			this.x0 = x0;
			this.t1 = t1;
			this.simulator = simulator;
		}

		public PDMPTrajectory simulate() {
			double[] x1 = new double[x0.length];
			PDMPTrajectory cm = new PDMPTrajectory();
			simulator.addStepHandler(cm);
			simulator.addReactionHandler(cm);
			simulator.simulate(model, t0, x0, t1, x1);
			simulator.clearReactionHandlers();
			simulator.clearStepHandlers();
			return cm;
		}

		@Override
		public PDMPTrajectory call() {
			return simulate();
		}

	}

	private class FiniteSimulationWorker implements Callable<FinitePDMPTrajectory> {

		private FinitePDMPTrajectory mt;
		private PDMPModel model;
		private double t0;
		private double[] x0;
		private double t1;
		private PDMPSimulator simulator;

		public FiniteSimulationWorker(PDMPSimulator simulator,
				FinitePDMPTrajectory mt, PDMPModel model, double t0,
				double[] x0, double t1) {
			this.mt = mt;
			this.model = model;
			this.t0 = t0;
			this.x0 = x0;
			this.t1 = t1;
			this.simulator = simulator;
		}

		public FinitePDMPTrajectory simulate() {
			double[] x1 = new double[x0.length];
			simulator.addStepHandler(mt);
			simulator.addReactionHandler(mt);
			simulator.simulate(model, t0, x0, t1, x1);
			simulator.clearReactionHandlers();
			simulator.clearStepHandlers();
			return mt;
		}

		@Override
		public FinitePDMPTrajectory call() {
			return simulate();
		}

	}

	private class TrajectoryDistributionSimulator extends SimulationWorker implements Runnable {

		private double[] tSeries;
		private SynchronizedSummaryStatistics[][] xStatistics;

		public TrajectoryDistributionSimulator(PDMPSimulator simulator,
				PDMPModel model, SynchronizedSummaryStatistics[][] xStatistics,
				double[] tSeries, double[] x0) {
			super(simulator, model, tSeries[0], x0, tSeries[tSeries.length - 1]);
			this.xStatistics = xStatistics;
			this.tSeries = tSeries;
		}

		@Override
		public void run() {
			PDMPTrajectory cm = simulate();
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

		public FiniteTrajectoryDistributionSimulator(
				PDMPSimulator simulator, FinitePDMPTrajectory mt,
				PDMPModel model, SynchronizedSummaryStatistics[][] xStatistics,
				double t0, double[] x0, double t1) {
			super(simulator, mt, model, t0, x0, t1);
			this.xStatistics = xStatistics;
		}

		@Override
		public void run() {
			FinitePDMPTrajectory mt = simulate();
			double[][] xSeries = mt.getxSeries();
			for (int s = 0; s < xStatistics.length; s++)
				for (int i = 0; i < xStatistics[s].length; i++)
					xStatistics[s][i].addValue(xSeries[s][i]);
		}

	}

	public final static int DEFAULT_NUMBER_OF_THREADS = 4;

	private PDMPSimulatorFactory simulatorFactory;
	private IntegratorFactory integratorFactory;
	private RandomDataGeneratorFactory rdgFactory;
	private ExecutorService executor;

	public PDMPSimulatorController() {
		this(DEFAULT_NUMBER_OF_THREADS);
	}

	public PDMPSimulatorController(int numOfThreads) {
		this(Executors.newFixedThreadPool(numOfThreads));
	}

	public PDMPSimulatorController(ExecutorService executor) {
		simulatorFactory = new DefaultPDMPSimulatorFactory();
		integratorFactory = new DefaultIntegratorFactory();
		rdgFactory = new DefaultRandomDataGeneratorFactory();
		this.executor = executor;
	}

	public void setExecutorService(ExecutorService executor) {
		this.executor = executor;
	}

	public void setSimulatorFactory(PDMPSimulatorFactory simulatorfactory) {
		this.simulatorFactory = simulatorfactory;
	}

	public void setIntegratorFactory(IntegratorFactory integratorFactory) {
		this.integratorFactory = integratorFactory;
	}

	public void setRandomDataGeneratorFactory(RandomDataGeneratorFactory rdgFactory) {
		this.rdgFactory = rdgFactory;
	}

	private PDMPSimulator createSimulator() {
		return simulatorFactory.createSimulator(createIntegrator(), createRandomDataGenerator());
	}

	private AbstractIntegrator createIntegrator() {
		return integratorFactory.createIntegrator();
	}

	private RandomDataGenerator createRandomDataGenerator() {
		return rdgFactory.createRandomDataGenerator();
	}

	public PDMPTrajectory simulateTrajectory(PDMPModel model, double t0, double[] x0, double t1) {
		PDMPSimulator simulator = createSimulator();
		SimulationWorker sw = new SimulationWorker(simulator, model, t0, x0, t1);
		return sw.simulate();
	}

	public double[][] simulateTrajectory(PDMPModel model, double[] tSeries, double[] x0) {
		FinitePDMPTrajectory mt = new FinitePDMPTrajectory(tSeries);
		simulateTrajectory(model, mt, tSeries[0], x0, tSeries[tSeries.length - 1]);
		return mt.getxSeries();
	}

	public void simulateTrajectory(PDMPModel model, FinitePDMPTrajectory mt, double t0, double[] x0, double t1) {
		PDMPSimulator simulator = createSimulator();
		FiniteSimulationWorker sw = new FiniteSimulationWorker(simulator, mt, model, t0, x0, t1);
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
			PDMPSimulator simulator = createSimulator();
			// PDMPModel has an internal state so we create a copy for each run
			PDMPModel model = modelFactory.createModel();
			Runnable r = new TrajectoryDistributionSimulator(simulator, model, xSeriesStatistics, tSeries, x0);
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
			PDMPSimulator simulator = createSimulator();
			// PDMPModel has an internal state so we create a copy for each run
			PDMPModel model = modelFactory.createModel();
			FinitePDMPTrajectory mt = mtFactory.createModelTrajectory();
			if (xSeriesStatistics == null) {
				xSeriesStatistics = new SynchronizedSummaryStatistics[x0.length][mt.gettSeries().length];
				for (int s = 0; s < xSeriesStatistics.length; s++)
					for (int i = 0; i < xSeriesStatistics[s].length; i++)
						xSeriesStatistics[s][i] = new SynchronizedSummaryStatistics();
			}
			Runnable r = new FiniteTrajectoryDistributionSimulator(simulator,
					mt, model, xSeriesStatistics, t0, x0.clone(), t1);
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
