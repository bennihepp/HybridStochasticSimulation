package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Semaphore;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.stat.descriptive.SynchronizedSummaryStatistics;

import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.providers.ObjProvider;
import ch.ethz.khammash.hybridstochasticsimulation.simulators.Simulator;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.DummyTrajectoryMapper;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteStatisticalSummaryTrajectory;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectory;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectoryMapper;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.VectorFiniteDistributionTrajectoryBuilder;

import com.google.common.base.Optional;
import com.google.common.util.concurrent.FutureCallback;
import com.google.common.util.concurrent.Futures;
import com.google.common.util.concurrent.ListenableFuture;
import com.google.common.util.concurrent.ListeningExecutorService;
import com.google.common.util.concurrent.MoreExecutors;

public abstract class AbstractSimulationController<T extends ReactionNetworkModel>
		implements SimulationController<T> {

	public final static int DEFAULT_NUMBER_OF_THREADS = 4;

	private static final int MAX_QUEUED_JOBS = 100;

	private ListeningExecutorService executor;
	private Optional<ObjProvider<? extends Simulator<T>>> simulatorProviderOptional;

	public AbstractSimulationController() {
		this(DEFAULT_NUMBER_OF_THREADS);
	}

	public AbstractSimulationController(int numOfThreads) {
		this(Executors.newFixedThreadPool(numOfThreads));
	}

	public AbstractSimulationController(ExecutorService executor) {
		this.executor = MoreExecutors.listeningDecorator(executor);
		simulatorProviderOptional = Optional.absent();
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

	protected Optional<ObjProvider<? extends Simulator<T>>> getSimulatorFactoryOptional() {
		return simulatorProviderOptional;
	}

	@Override
	final public void setSimulatorProvider(ObjProvider<? extends Simulator<T>> simulatorProvider) {
		simulatorProviderOptional = Optional.<ObjProvider<? extends Simulator<T>>>of(simulatorProvider);
	}

//	@Override
//	final public void setTrajectoryRecorderFactory(
//			TrajectoryRecorderFactory<E> trajectoryRecorderFactory) {
//		trajectoryRecorderFactoryOptional = Optional.of(trajectoryRecorderFactory);
//	}

	protected Simulator<T> createSimulator() {
        // We want to have different random number sequences for each run
		if (simulatorProviderOptional.isPresent())
			return simulatorProviderOptional.get().get();
		else
			throw new UnsupportedOperationException("The SimulatorProvider has not been set yet");
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
	public List<TrajectoryRecorder> simulateTrajectories(
			int runs, ObjProvider<? extends T> modelProvider, ObjProvider<? extends TrajectoryRecorder> trProvider,
			double t0, double[] x0, double t1) {
		checkArgument(runs > 0, "Expected runs > 0");

		final List<TrajectoryRecorder> trList = Collections.synchronizedList(new LinkedList<TrajectoryRecorder>());

        final Semaphore jobQueueCounter = new Semaphore(MAX_QUEUED_JOBS);
//        final AtomicInteger queuedJobs = new AtomicInteger(0);

        for (int k=0; k < runs; k++) {
        	jobQueueCounter.acquireUninterruptibly();
//        	makeSubmission(queuedJobs);
    		T model = modelProvider.get();
    		TrajectoryRecorder tr = trProvider.get();
    		Callable<TrajectoryRecorder> sw = createSimulationWorker(model, tr, t0, x0, t1);
    		ListenableFuture<TrajectoryRecorder> completion = executor.submit(sw);
    		Futures.addCallback(completion, new FutureCallback<TrajectoryRecorder>() {

				@Override
				public void onFailure(Throwable t) {
					jobQueueCounter.release();
//					releaseSubmission(queuedJobs);
					throw new RuntimeException(t);
				}

				@Override
				public void onSuccess(TrajectoryRecorder tr) {
					jobQueueCounter.release();
//					releaseSubmission(queuedJobs);
					trList.add(tr);
				}

    		});
        }
        executor.shutdown();
        do {
        	try {
        		executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
	        } catch (InterruptedException e) { }
        } while (!executor.isTerminated());
        return trList;
	}

	@Override
	public FiniteStatisticalSummaryTrajectory simulateTrajectoryDistribution(
			int runs, ObjProvider<? extends T> modelProvider, ObjProvider<? extends FiniteTrajectoryRecorder> trProvider,
			double t0, double[] x0, double t1) {
		DummyTrajectoryMapper mapper = new DummyTrajectoryMapper();
		return simulateTrajectoryDistribution(runs, modelProvider, trProvider, mapper, t0, x0, t1);
	}

	@Override
	public FiniteStatisticalSummaryTrajectory simulateTrajectoryDistribution(
			int runs, ObjProvider<? extends T> modelProvider, ObjProvider<? extends FiniteTrajectoryRecorder> trProvider,
			final FiniteTrajectoryMapper mapper, double t0, double[] x0, double t1) {
		checkArgument(runs > 0, "Expected runs > 0");
//		checkArgument(tSeries.length >= 2, "Expected tSeries.length >= 2");

//		FiniteTrajectoryRecorder dummyTr = trProvider.get();
//		dummyTr.beginRecording(t0, x0, t1);
//		int numOfStates = dummyTr.getNumberOfStates();
//		int numOfTimePoints = dummyTr.getNumberOfTimePoints();
//		double[] tSeries = dummyTr.gettSeries();
//		final SynchronizedSummaryStatistics[][] xSeriesStatistics = new SynchronizedSummaryStatistics[numOfStates][numOfTimePoints];
//        for (int s=0; s < numOfStates; s++)
//        	for (int i=0; i < dummyTr.getNumberOfTimePoints(); i++)
//				xSeriesStatistics[s][i] = new SynchronizedSummaryStatistics();
		final VectorFiniteDistributionTrajectoryBuilder trajectoryDistributionBuilder = new VectorFiniteDistributionTrajectoryBuilder();

        final Semaphore jobQueueCounter = new Semaphore(MAX_QUEUED_JOBS);
//        final AtomicInteger queuedJobs = new AtomicInteger(0);

        for (int k=0; k < runs; k++) {
        	jobQueueCounter.acquireUninterruptibly();
//        	makeSubmission(queuedJobs);
    		T model = modelProvider.get();
    		FiniteTrajectoryRecorder tr = trProvider.get();
    		Callable<TrajectoryRecorder> sw = createSimulationWorker(model, tr, t0, x0, t1);
    		ListenableFuture<TrajectoryRecorder> completion = executor.submit(sw);
    		Futures.addCallback(completion, new FutureCallback<TrajectoryRecorder>() {

				@Override
				public void onFailure(Throwable t) {
					jobQueueCounter.release();
//					releaseSubmission(queuedJobs);
					throw new RuntimeException(t);
				}

				@Override
				public void onSuccess(TrajectoryRecorder tr) {
					jobQueueCounter.release();
//					releaseSubmission(queuedJobs);
					FiniteTrajectoryRecorder ftr = (FiniteTrajectoryRecorder)tr;
					FiniteTrajectory mappedTr = mapper.map(ftr);
					trajectoryDistributionBuilder.addTrajectory(mappedTr);
//					double[][] xSeries = ftr.getxSeries();
//					for (int s = 0; s < ftr.getNumberOfStates(); s++) {
//						for (int i = 0; i < ftr.getNumberOfTimePoints(); i++) {
//							xSeriesStatistics[s][i].addValue(xSeries[s][i]);
//						}
//					}
				}

    		});
        }
        executor.shutdown();
        do {
        	try {
        		executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
	        } catch (InterruptedException e) { }
        } while (!executor.isTerminated());
        return trajectoryDistributionBuilder.getDistributionTrajectory();
//        return VectorFiniteStatisticalSummaryPlotData.createFromStatisticalSummary(tSeries, xSeriesStatistics);
	}

//	private void makeSubmission(AtomicInteger queuedJobs) {
//		synchronized (queuedJobs) {
//        	while (queuedJobs.get() >= MAX_QUEUED_JOBS) {
//				try {
//					queuedJobs.wait();
//				} catch (InterruptedException e) {}
//        	}
//    		queuedJobs.getAndIncrement();
//		}
//	}

//	private void releaseSubmission(AtomicInteger queuedJobs) {
//		synchronized(queuedJobs) {
//			queuedJobs.getAndDecrement();
//			queuedJobs.notifyAll();
//		}
//	}

}
