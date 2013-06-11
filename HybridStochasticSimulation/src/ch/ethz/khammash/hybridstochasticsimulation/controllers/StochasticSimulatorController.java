package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import static com.google.common.base.Preconditions.checkArgument;

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

import ch.ethz.khammash.hybridstochasticsimulation.models.UnaryBinaryStochasticModel;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.StochasticSimulator;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteStochasticTrajectory;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.StochasticTrajectory;


public class StochasticSimulatorController {

	private class SimulationWorker implements Callable<StochasticTrajectory> {

		protected UnaryBinaryStochasticModel model;
		protected double t0;
		protected double[] x0;
		protected double t1;
		protected StochasticSimulator simulator;

		public SimulationWorker(RandomDataGenerator rng, UnaryBinaryStochasticModel model, double t0, double[] x0, double t1) {
			this.model = model;
			this.t0 = t0;
			this.x0 = x0;
			this.t1 = t1;
			simulator = new StochasticSimulator(rng);
		}

		public StochasticTrajectory simulate() {
			double[] x1 = new double[x0.length];
			StochasticTrajectory smt = new StochasticTrajectory();
			simulator.addReactionHandler(smt);
			simulator.simulate(model, t0, x0, t1, x1);
			simulator.clearReactionHandlers();
			return smt;
		}

		@Override
		public StochasticTrajectory call() {
			return simulate();
		}

	}

	private class FiniteSimulationWorker implements Callable<FiniteStochasticTrajectory> {

		protected UnaryBinaryStochasticModel model;
		protected double[] tSeries;
		protected double[] x0;
		protected StochasticSimulator simulator;

		public FiniteSimulationWorker(RandomDataGenerator rng, UnaryBinaryStochasticModel model, double[] tSeries, double[] x0) {
			this.model = model;
			this.tSeries = tSeries;
			this.x0 = x0;
			simulator = new StochasticSimulator(rng);
		}

		public FiniteStochasticTrajectory simulate() {
			double[] x1 = new double[x0.length];
			FiniteStochasticTrajectory smt = new FiniteStochasticTrajectory(tSeries);
			simulator.addReactionHandler(smt);
			simulator.simulate(model, tSeries[0], x0, tSeries[tSeries.length - 1], x1);
			simulator.clearReactionHandlers();
			return smt;
		}

		@Override
		public FiniteStochasticTrajectory call() {
			return simulate();
		}

	}

	private class TrajectoryDistributionSimulator extends FiniteSimulationWorker implements Runnable {

		protected double[] tSeries;
		protected SynchronizedSummaryStatistics[][] xStatistics;

		TrajectoryDistributionSimulator(RandomDataGenerator rng, UnaryBinaryStochasticModel model,
				SynchronizedSummaryStatistics[][] xStatistics, double[] tSeries, double[] x0) {
			super(rng, model, tSeries, x0);
			this.xStatistics = xStatistics;
			this.tSeries = tSeries;

		}

		@Override
		public void run() {
			FiniteStochasticTrajectory smt = simulate();
			for (int i = 0; i < tSeries.length; i++) {
				double[] x = smt.getInterpolatedState(tSeries[i]);
				for (int s = 0; s < x.length; s++) {
					try {
						xStatistics[i][s].addValue(x[s]);
					} catch (Exception e) {
						e.printStackTrace();
					}
				}
			}
		}
	}

    
    public final static int DEFAULT_NUMBER_OF_THREADS = 4;

    private ExecutorService executor;
    private RandomGenerator rng;

    public StochasticSimulatorController() {
        this(DEFAULT_NUMBER_OF_THREADS);
    }

	public StochasticSimulatorController(int numOfThreads) {
        executor = Executors.newFixedThreadPool(numOfThreads);
    }

	public StochasticSimulatorController(ExecutorService executor) {
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
		checkArgument(numOfGenerators > 0, "Expected numOfGenerators > 0");
    	RandomDataGenerator[] result = new RandomDataGenerator[numOfGenerators];
    	for (int i=0; i < result.length; i++) {
    		result[i] = createRandomDataGenerator();
    	}
    	return result;
    }

	public FiniteStochasticTrajectory simulateTrajectory(UnaryBinaryStochasticModel model, double[] tSeries, double[] x0) {
		checkArgument(tSeries.length >= 2, "Expected tSeries.length >= 2");
		RandomDataGenerator rdg = createRandomDataGenerator();
    	FiniteSimulationWorker sw = new FiniteSimulationWorker(rdg, model, tSeries, x0);
    	FiniteStochasticTrajectory mt = sw.simulate();
		return mt;
//		StochasticModelTrajectory mt = simulateFixedTrajectory(tSeries[0], x0, tSeries[tSeries.length - 1]);
//    	double[][] xSeries = new double[tSeries.length][x0.length];
//		for (int i = 0; i < tSeries.length; i++) {
//			double[] x = mt.getInterpolatedState(tSeries[i]);
//			xSeries[i] = x;
//		}
//        return xSeries;
	}

	public StochasticTrajectory simulateTrajectory(UnaryBinaryStochasticModel model, double t0, double[] x0, double t1) {
		RandomDataGenerator rdg = createRandomDataGenerator();
    	SimulationWorker sw = new SimulationWorker(rdg, model, t0, x0, t1);
		StochasticTrajectory mt = sw.simulate();
		return mt;
	}

	public StatisticalSummary[][] computeTrajectoryDistribution(UnaryBinaryStochasticModel model, int runs, double[] tSeries, double[] x0)
			throws InterruptedException, CancellationException, ExecutionException {
		checkArgument(runs > 0, "Expected runs > 0");
		checkArgument(tSeries.length >= 2, "Expected tSeries.length >= 2");
		SynchronizedSummaryStatistics[][] xSeriesStatistics = new SynchronizedSummaryStatistics[tSeries.length][x0.length];
        for (int i=0; i < tSeries.length; i++)
	        for (int s=0; s < x0.length; s++)
				xSeriesStatistics[i][s] = new SynchronizedSummaryStatistics();
        // We want to have different random number sequences for each run
        RandomDataGenerator[] rdgs = createRandomDataGenerators(runs);
        for (int k=0; k < runs; k++) {
			Runnable r = new TrajectoryDistributionSimulator(rdgs[k], model, xSeriesStatistics, tSeries, x0);
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
