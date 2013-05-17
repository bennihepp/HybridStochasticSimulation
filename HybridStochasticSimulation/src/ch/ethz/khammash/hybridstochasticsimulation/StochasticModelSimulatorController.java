package ch.ethz.khammash.hybridstochasticsimulation;

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

/**
 * Created by
 * User: bhepp
 * Date: 5/13/13
 * Time: 10:51 AM
 */

public class StochasticModelSimulatorController {

	private class SimulationWorker implements Callable<StochasticModelTrajectory> {

		protected StochasticModel model;
		protected double t0;
		protected double[] x0;
		protected double t1;
		protected StochasticModelSimulator simulator;

		public SimulationWorker(RandomDataGenerator rng, StochasticModel model, double t0, double[] x0, double t1) {
			this.model = model;
			this.t0 = t0;
			this.x0 = x0;
			this.t1 = t1;
			simulator = new StochasticModelSimulator(rng);
		}

		public StochasticModelTrajectory simulate() {
			double[] x1 = new double[x0.length];
			StochasticModelTrajectory smt = new StochasticModelTrajectory();
			simulator.addReactionHandler(smt);
			simulator.simulate(model, t0, x0, t1, x1);
			simulator.clearReactionHandlers();
			return smt;
		}

		@Override
		public StochasticModelTrajectory call() {
			return simulate();
		}

	}

	private class TrajectoryDistributionSimulator extends SimulationWorker implements Runnable {

		protected double[] tSeries;
		protected SynchronizedSummaryStatistics[][] xStatistics;

		TrajectoryDistributionSimulator(RandomDataGenerator rng, StochasticModel model,
				SynchronizedSummaryStatistics[][] xStatistics, double[] tSeries, double[] x0) {
			super(rng, model, tSeries[0], x0, tSeries[tSeries.length - 1]);
			this.xStatistics = xStatistics;
			this.tSeries = tSeries;

		}

		@Override
		public void run() {
			StochasticModelTrajectory smt = simulate();
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
    private StochasticModel model;
    private RandomGenerator rng;

    public StochasticModelSimulatorController(StochasticModel model) {
        this(model, DEFAULT_NUMBER_OF_THREADS);
    }

	public StochasticModelSimulatorController(StochasticModel model, int numOfThreads) {
        this.model = model;
        executor = Executors.newFixedThreadPool(numOfThreads);
    }

	public StochasticModelSimulatorController(StochasticModel mode, ExecutorService executor) {
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
    	for (int i=0; i < result.length; i++) {
    		result[i] = createRandomDataGenerator();
    	}
    	return result;
    }

	public double[][] simulateTrajectory(double[] tSeries, double[] x0) {
		StochasticModelTrajectory mt = simulateTrajectory(tSeries[0], x0, tSeries[tSeries.length - 1]);
    	double[][] xSeries = new double[tSeries.length][x0.length];
		for (int i = 0; i < tSeries.length; i++) {
			double[] x = mt.getInterpolatedState(tSeries[i]);
			xSeries[i] = x;
		}
        return xSeries;
	}

	public StochasticModelTrajectory simulateTrajectory(double t0, double[] x0, double t1) {
		RandomDataGenerator rdg = createRandomDataGenerator();
    	SimulationWorker sw = new SimulationWorker(rdg, model, t0, x0, t1);
		StochasticModelTrajectory mt = sw.simulate();
		return mt;
	}

	public StatisticalSummary[][] computeTrajectoryDistribution(int runs, double[] tSeries, double[] x0)
			throws InterruptedException, CancellationException, ExecutionException {
		SynchronizedSummaryStatistics[][] xSeriesStatistics = new SynchronizedSummaryStatistics[tSeries.length][x0.length];
        for (int i=0; i < tSeries.length; i++)
	        for (int s=0; s < x0.length; s++)
				xSeriesStatistics[i][s] = new SynchronizedSummaryStatistics();
        // We want to have different random number sequences for each run
        RandomDataGenerator[] rdgs = createRandomDataGenerators(runs);
        final long startTime = System.currentTimeMillis();
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
        final long endTime = System.currentTimeMillis();
        System.out.println("Execution time: " + (endTime - startTime) );
        return xSeriesStatistics;
    }

}
