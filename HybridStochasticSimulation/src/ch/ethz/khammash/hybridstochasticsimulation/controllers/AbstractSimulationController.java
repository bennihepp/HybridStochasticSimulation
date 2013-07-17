package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.concurrent.Callable;
import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.stat.descriptive.SynchronizedSummaryStatistics;

import ch.ethz.khammash.hybridstochasticsimulation.factories.DefaultRandomDataGeneratorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.FiniteTrajectoryRecorderFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.ModelFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.RandomDataGeneratorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.factories.SimulatorFactory;
import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteStatisticalSummaryTrajectory;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.VectorFiniteStatisticalSummaryPlotData;

import com.google.common.base.Optional;
import com.google.common.util.concurrent.FutureCallback;
import com.google.common.util.concurrent.Futures;
import com.google.common.util.concurrent.ListenableFuture;
import com.google.common.util.concurrent.ListeningExecutorService;
import com.google.common.util.concurrent.MoreExecutors;

public abstract class AbstractSimulationController<T extends ReactionNetworkModel>
		implements SimulationController<T> {

	public final static int DEFAULT_NUMBER_OF_THREADS = 4;

	private ListeningExecutorService executor;
	private Optional<SimulatorFactory<Simulator<T>>> simulatorFactoryOptional;
//	private Optional<TrajectoryRecorderFactory<E>> trajectoryRecorderFactoryOptional;
	private Optional<RandomDataGeneratorFactory> rdgFactoryOptional;

	public AbstractSimulationController() {
		this(DEFAULT_NUMBER_OF_THREADS);
	}

	public AbstractSimulationController(int numOfThreads) {
		this(Executors.newFixedThreadPool(numOfThreads));
	}

	public AbstractSimulationController(ExecutorService executor) {
		this.executor = MoreExecutors.listeningDecorator(executor);
		simulatorFactoryOptional = Optional.absent();
//		trajectoryRecorderFactoryOptional = Optional.absent();
		rdgFactoryOptional = Optional.<RandomDataGeneratorFactory>of(new DefaultRandomDataGeneratorFactory());
	}

	protected SimulationWorker createSimulationWorker(T model, TrajectoryRecorder tr, double t0, double[] x0, double t1) {
		Simulator<T> simulator = createSimulator();
		return new DefaultSimulationWorker<T>(simulator, model, tr, t0, x0, t1);
	}

	protected FiniteSimulationWorker createFiniteSimulationWorker(T model, FiniteTrajectoryRecorder tr, double t0, double[] x0, double t1) {
		Simulator<T> simulator = createSimulator();
		return new DefaultFiniteSimulationWorker<T>(simulator, model, tr, t0, x0, t1);
	}

	protected DistributionSimulationWorker createDistributionSimulationWorker(
			T model, FiniteTrajectoryRecorder tr, double[] tSeries, SynchronizedSummaryStatistics[][] xSeriesStatistics,
			double t0, double[] x0, double t1) {
//		checkArgument(tSeries.length >= 2, "Expected tSeries.length >= 2");
        FiniteSimulationWorker worker = createFiniteSimulationWorker(model, tr, t0, x0, t1);
		return new DefaultDistributionSimulationWorker(worker, tSeries, xSeriesStatistics);
	}

	@Override
	final public void setExecutorService(ExecutorService executor) {
		this.executor = MoreExecutors.listeningDecorator(executor);
	}

	protected Optional<SimulatorFactory<Simulator<T>>> getSimulatorFactoryOptional() {
		return simulatorFactoryOptional;
	}

	@Override
	final public void setSimulatorFactory(SimulatorFactory<Simulator<T>> simulatorFactory) {
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

	protected Simulator<T> createSimulator() {
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
	public void simulateTrajectory(T model, TrajectoryRecorder tr, double t0, double[] x0, double t1) {
		SimulationWorker sw = createSimulationWorker(model, tr, t0, x0, t1);
		sw.simulate();
	}

	@Override
	public FiniteStatisticalSummaryTrajectory simulateTrajectoryDistribution(
			int runs, ModelFactory<T> modelFactory, FiniteTrajectoryRecorderFactory trFactory,
			double t0, double[] x0, double t1)
			throws InterruptedException, CancellationException, ExecutionException {
		checkArgument(runs > 0, "Expected runs > 0");
//		checkArgument(tSeries.length >= 2, "Expected tSeries.length >= 2");

		FiniteTrajectoryRecorder dummyTr = trFactory.createTrajectoryRecorder();
		dummyTr.beginRecording(t0, x0, t1);
		int numOfStates = dummyTr.getNumberOfStates();
		int numOfTimePoints = dummyTr.getNumberOfTimePoints();
		double[] tSeries = dummyTr.gettSeries();
		final SynchronizedSummaryStatistics[][] xSeriesStatistics = new SynchronizedSummaryStatistics[numOfStates][numOfTimePoints];
        for (int s=0; s < numOfStates; s++)
        	for (int i=0; i < dummyTr.getNumberOfTimePoints(); i++)
				xSeriesStatistics[s][i] = new SynchronizedSummaryStatistics();

        for (int k=0; k < runs; k++) {
    		T model = modelFactory.createModel();
    		FiniteTrajectoryRecorder tr = trFactory.createTrajectoryRecorder();
    		Callable<TrajectoryRecorder> sw = createSimulationWorker(model, tr, t0, x0, t1);
    		ListenableFuture<TrajectoryRecorder> completion = executor.submit(sw);
    		Futures.addCallback(completion, new FutureCallback<TrajectoryRecorder>() {

				@Override
				public void onFailure(Throwable t) {
					throw new RuntimeException(t);
				}

				@Override
				public void onSuccess(TrajectoryRecorder tr) {
					FiniteTrajectoryRecorder ftr = (FiniteTrajectoryRecorder)tr;
					double[][] xSeries = ftr.getxSeries();
					for (int s = 0; s < ftr.getNumberOfStates(); s++) {
						for (int i = 0; i < ftr.getNumberOfTimePoints(); i++) {
							xSeriesStatistics[s][i].addValue(xSeries[s][i]);
						}
					}
				}

    		});
        }
        executor.shutdown();
        do {
        	try {
        		executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
	        } catch (InterruptedException e) { }
        } while (!executor.isTerminated());
        return VectorFiniteStatisticalSummaryPlotData.createFromStatisticalSummary(tSeries, xSeriesStatistics);
	}

}
