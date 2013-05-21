package ch.ethz.khammash.hybridstochasticsimulation.simulators;

import java.util.concurrent.Callable;
import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.StatisticalSummary;
import org.apache.commons.math3.stat.descriptive.SynchronizedSummaryStatistics;

import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;

/**
 * Created by
 * User: bhepp
 * Date: 5/13/13
 * Time: 10:51 AM
 */

public class PDMPModelSimulatorController {

	private class SimulationWorker implements Callable<PDMPModelTrajectory> {

		protected PDMPModel model;
		protected double t0;
		protected double[] x0;
		protected double t1;
		protected PDMPModelSimulator simulator;

		public SimulationWorker(RandomDataGenerator rng, PDMPModel model, double t0, double[] x0, double t1) {
			this.model = model;
			this.t0 = t0;
			this.x0 = x0;
			this.t1 = t1;
			simulator = new PDMPModelSimulator(rng);
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

	private class TrajectoryDistributionSimulator extends SimulationWorker implements Runnable {
		protected double[] tSeries;
		protected SynchronizedSummaryStatistics[][] xStatistics;

		TrajectoryDistributionSimulator(RandomDataGenerator rng, PDMPModel model, SynchronizedSummaryStatistics[][] xStatistics,
				double[] tSeries, double[] x0) {
			super(rng, model, tSeries[0], x0, tSeries[tSeries.length - 1]);
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
					// zDistribution[i][s] += z[s];
					xStatistics[i][s].addValue(x[s]);
				}
			}
		}
	}

	public final static int DEFAULT_NUMBER_OF_THREADS = 4;

	private ExecutorService executor;
	private PDMPModel model;
	private RandomGenerator rng;

	public PDMPModelSimulatorController(PDMPModel model) {
		this(model, DEFAULT_NUMBER_OF_THREADS);
	}

	public PDMPModelSimulatorController(PDMPModel model, int numOfThreads) {
		this.model = model;
		executor = Executors.newFixedThreadPool(numOfThreads);
	}

	public PDMPModelSimulatorController(PDMPModel mode, ExecutorService executor) {
		this.executor = executor;
	}

	public void setRandomGenerator(RandomGenerator rng) {
		this.rng = rng;
	}

	private RandomDataGenerator createRandomDataGenerator() {
		if (rng == null)
			rng = new MersenneTwister();
		RandomDataGenerator rdg = new RandomDataGenerator();
		rdg.reSeed(rng.nextLong());
		return rdg;
	}

	private RandomDataGenerator[] createRandomDataGenerators(int numOfGenerators) {
		RandomDataGenerator[] result = new RandomDataGenerator[numOfGenerators];
		for (int i = 0; i < result.length; i++) {
			result[i] = createRandomDataGenerator();
		}
		return result;
	}

	public double[][] simulateTrajectory(double[] tSeries, double[] x0) {
		RandomDataGenerator rdg = createRandomDataGenerator();
		double[][] xSeries = new double[tSeries.length][x0.length];
		SimulationWorker sw = new SimulationWorker(rdg, model, tSeries[0], x0, tSeries[tSeries.length - 1]);

		PDMPModelTrajectory mt = sw.simulate();
		for (int i = 0; i < tSeries.length; i++) {
			double[] x = mt.getInterpolatedState(tSeries[i]);
			xSeries[i] = x;
		}
		return xSeries;
	}

	public StatisticalSummary[][] simulateTrajectoryDistribution(int runs, double[] tSeries, double[] x0)
			throws InterruptedException, CancellationException, ExecutionException {
		SynchronizedSummaryStatistics[][] xSeriesStatistics = new SynchronizedSummaryStatistics[tSeries.length][x0.length];
		for (int i = 0; i < tSeries.length; i++)
			for (int s = 0; s < x0.length; s++)
				xSeriesStatistics[i][s] = new SynchronizedSummaryStatistics();
		// PDMPModel has an internal state so we create a copy for each run
		PDMPModel[] models = new PDMPModel[runs];
		for (int k = 0; k < runs; k++)
			models[k] = new PDMPModel(model);
		// We also want to have different random number sequences for each run
		RandomDataGenerator[] rdgs = createRandomDataGenerators(runs);
		final long startTime = System.currentTimeMillis();
		for (int k = 0; k < runs; k++) {
			Runnable r = new TrajectoryDistributionSimulator(rdgs[k], models[k], xSeriesStatistics, tSeries, x0);
			executor.execute(r);
		}
		executor.shutdown();
		try {
			executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
		} catch (InterruptedException e) {
			executor.shutdownNow();
			Thread.currentThread().interrupt();
		}
		final long endTime = System.currentTimeMillis();
		System.out.println("Execution time: " + (endTime - startTime));
		return xSeriesStatistics;
	}

}
