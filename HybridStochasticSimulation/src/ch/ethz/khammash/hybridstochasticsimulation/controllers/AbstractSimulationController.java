package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.stat.descriptive.StatisticalSummary;
import org.apache.commons.math3.stat.descriptive.SynchronizedSummaryStatistics;

import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;

import com.google.common.base.Optional;

public abstract class AbstractSimulationController<T extends ReactionNetworkModel, E extends TrajectoryRecorder<T>>
		implements SimulationController<T, E> {

	public final static int DEFAULT_NUMBER_OF_THREADS = 4;

	private ExecutorService executor;
	private Optional<SimulatorFactory<Simulator<T, E>>> simulatorFactoryOptional;
//	private Optional<TrajectoryRecorderFactory<E>> trajectoryRecorderFactoryOptional;
	private Optional<RandomDataGeneratorFactory> rdgFactoryOptional;

	public AbstractSimulationController() {
		this(DEFAULT_NUMBER_OF_THREADS);
	}

	public AbstractSimulationController(int numOfThreads) {
		this(Executors.newFixedThreadPool(numOfThreads));
	}

	public AbstractSimulationController(ExecutorService executor) {
		this.executor = executor;
		simulatorFactoryOptional = Optional.absent();
//		trajectoryRecorderFactoryOptional = Optional.absent();
		rdgFactoryOptional = Optional.<RandomDataGeneratorFactory>of(new DefaultRandomDataGeneratorFactory());
	}

	protected SimulationWorker<T, E> createSimulationWorker(T model, E tr, double t0, double[] x0, double t1) {
		Simulator<T, E> simulator = createSimulator();
		return new DefaultSimulationWorker<T, E>(simulator, model, tr, t0, x0, t1);
	}

	protected DistributionSimulationWorker createDistributionSimulationWorker(
			T model, E tr, SynchronizedSummaryStatistics[][] xSeriesStatistics, double[] tSeries, double[] x0) {
		checkArgument(tSeries.length >= 2, "Expected tSeries.length >= 2");
		double t0 = tSeries[0];
		double t1 = tSeries[tSeries.length - 1];
        SimulationWorker<T, E> worker = createSimulationWorker(model, tr, t0, x0, t1);
		return new DefaultDistributionSimulationWorker<T, E>(worker, xSeriesStatistics, tSeries, x0);
	}

	@Override
	final public void setExecutorService(ExecutorService executor) {
		this.executor = executor;
	}

	protected Optional<SimulatorFactory<Simulator<T,E>>> getSimulatorFactoryOptional() {
		return simulatorFactoryOptional;
	}

	@Override
	final public void setSimulatorFactory(
			SimulatorFactory<Simulator<T, E>> simulatorFactory) {
		simulatorFactoryOptional = Optional.of(simulatorFactory);
	}

//	@Override
//	final public void setTrajectoryRecorderFactory(
//			TrajectoryRecorderFactory<E> trajectoryRecorderFactory) {
//		trajectoryRecorderFactoryOptional = Optional.of(trajectoryRecorderFactory);
//	}

	@Override
	final public void setRandomDataGeneratorFactory(RandomDataGeneratorFactory rdgFactory) {
		rdgFactoryOptional = Optional.of(rdgFactory);
	}

	protected Simulator<T, E> createSimulator() {
        // We want to have different random number sequences for each run
		if (simulatorFactoryOptional.isPresent())
			return simulatorFactoryOptional.get().createSimulator(createRandomDataGenerator());
		else
			throw new UnsupportedOperationException("The SimulatorFactory has not been set yet");
	}

//	private E createTrajectoryRecorder() {
//		if (trajectoryRecorderFactoryOptional.isPresent())
//			return trajectoryRecorderFactoryOptional.get().createTrajectoryRecorder();
//		else
//			throw new UnsupportedOperationException("The TrajectoryRecorderFactory has not been set yet");
//	}

//	private E createTrajectoryRecorder(double[] tSeries) {
//		if (trajectoryRecorderFactoryOptional.isPresent())
//			return trajectoryRecorderFactoryOptional.get().createTrajectoryRecorder(tSeries);
//		else
//			throw new UnsupportedOperationException("The TrajectoryRecorderFactory has not been set yet");
//	}

	protected RandomDataGenerator createRandomDataGenerator() {
		if (rdgFactoryOptional.isPresent())
			return rdgFactoryOptional.get().createRandomDataGenerator();
		else
			throw new UnsupportedOperationException("The RandomDataGeneratorFactory has not been set yet");
	}

//	@Override
//	public E simulateTrajectory(T model, double t0, double[] x0, double t1) {
//		E tr = createTrajectoryRecorder();
//		SimulationWorker<T, E> sw = createSimulationWorker(model, tr, t0, x0, t1);
//		return sw.simulate();
//	}

	@Override
	public void simulateTrajectory(T model, E tr, double t0, double[] x0, double t1) {
		SimulationWorker<T, E> sw = createSimulationWorker(model, tr, t0, x0, t1);
		sw.simulate();
	}

	@Override
	public StatisticalSummary[][] simulateTrajectoryDistribution(
			int runs, ModelFactory<T> modelFactory, FiniteTrajectoryRecorderFactory<E> trFactory,
			double[] tSeries, double[] x0)
			throws InterruptedException, CancellationException, ExecutionException {
		checkArgument(runs > 0, "Expected runs > 0");
		checkArgument(tSeries.length >= 2, "Expected tSeries.length >= 2");
		SynchronizedSummaryStatistics[][] xSeriesStatistics = new SynchronizedSummaryStatistics[x0.length][tSeries.length];
        for (int s=0; s < x0.length; s++)
        	for (int i=0; i < tSeries.length; i++)
				xSeriesStatistics[s][i] = new SynchronizedSummaryStatistics();
        for (int k=0; k < runs; k++) {
    		T model = modelFactory.createModel();
    		E tr = trFactory.createTrajectoryRecorder(tSeries);
            Runnable r = createDistributionSimulationWorker(model, tr, xSeriesStatistics, tSeries, x0);
            executor.execute(r);
        }
        executor.shutdown();
        do {
        	try {
        		executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
	        } catch (InterruptedException e) { }
        } while (!executor.isTerminated());
        return xSeriesStatistics;
	}

}
