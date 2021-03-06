package ch.ethz.bhepp.hybridstochasticsimulation.controllers;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Semaphore;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;

import org.apache.commons.math3.stat.descriptive.SynchronizedSummaryStatistics;
import org.apache.commons.math3.util.FastMath;

import ch.ethz.bhepp.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.bhepp.hybridstochasticsimulation.providers.ObjProvider;
import ch.ethz.bhepp.hybridstochasticsimulation.simulators.Simulator;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.DummyTrajectoryMapper;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FiniteStatisticalSummaryTrajectory;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FiniteTrajectory;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FiniteTrajectoryMapper;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.TrajectoryRecorder;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.VectorFiniteDistributionTrajectoryBuilder;

import com.google.common.util.concurrent.FutureCallback;
import com.google.common.util.concurrent.Futures;
import com.google.common.util.concurrent.ListenableFuture;
import com.google.common.util.concurrent.ListeningExecutorService;
import com.google.common.util.concurrent.MoreExecutors;

public abstract class AbstractSimulationController<T extends ReactionNetworkModel>
		implements SimulationController<T> {

	public static final int DEFAULT_NUMBER_OF_THREADS = 4;
	private static final int DEFAULT_MAX_QUEUED_JOBS = 100;

	private final ListeningExecutorService executor;
	private ObjProvider<? extends Simulator<T>> simulatorProvider;
	private int maxQueuedJobs = DEFAULT_MAX_QUEUED_JOBS;

	public AbstractSimulationController(ObjProvider<? extends Simulator<T>> simulatorProvider) {
		this(simulatorProvider, DEFAULT_NUMBER_OF_THREADS);
	}

	public AbstractSimulationController(ObjProvider<? extends Simulator<T>> simulatorProvider, int numOfThreads) {
		this(simulatorProvider, Executors.newFixedThreadPool(numOfThreads), FastMath.max(numOfThreads, DEFAULT_MAX_QUEUED_JOBS));
//		setMaxQueuedJobs(FastMath.max(numOfThreads, maxQueuedJobs));
	}

	public AbstractSimulationController(ObjProvider<? extends Simulator<T>> simulatorProvider, ExecutorService executor) {
		this(simulatorProvider, executor, DEFAULT_MAX_QUEUED_JOBS);
	}

	protected AbstractSimulationController(ObjProvider<? extends Simulator<T>> simulatorProvider, ExecutorService executor, int bound) {
		this.simulatorProvider = simulatorProvider;
		this.executor = MoreExecutors.listeningDecorator(BlockingThreadPoolExecutor.blockingDecorator(executor, bound));
//		this.executor = BlockingThreadPoolExecutor.blockingDecorator(MoreExecutors.listeningDecorator(executor), bound);
//		this.executor = MoreExecutors.listeningDecorator(executor);
//		simulatorProviderOptional = Optional.absent();
		setMaxQueuedJobs(bound);
	}

	public void setMaxQueuedJobs(int maxQueuedJobs) {
		checkArgument(maxQueuedJobs > 0, "Expected maxQuededJobs > 0, found maxQuededJobs == " + maxQueuedJobs);
		this.maxQueuedJobs = maxQueuedJobs;
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

//	@Override
//	final public void setExecutorService(ExecutorService executor) {
//		this.executor = MoreExecutors.listeningDecorator(executor);
//	}

//	protected Optional<ObjProvider<? extends Simulator<T>>> getSimulatorProviderOptional() {
//		return simulatorProviderOptional;
//	}

//	@Override
//	final public void setSimulatorProvider(ObjProvider<? extends Simulator<T>> simulatorProvider) {
//		simulatorProviderOptional = Optional.<ObjProvider<? extends Simulator<T>>>of(simulatorProvider);
//	}

//	@Override
//	final public void setTrajectoryRecorderFactory(
//			TrajectoryRecorderFactory<E> trajectoryRecorderFactory) {
//		trajectoryRecorderFactoryOptional = Optional.of(trajectoryRecorderFactory);
//	}

	protected Simulator<T> createSimulator() {
//		if (simulatorProviderOptional.isPresent())
//			return simulatorProviderOptional.get().get();
//		else
//			throw new UnsupportedOperationException("The SimulatorProvider has not been set yet");
		return simulatorProvider.get();
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
			double t0, double[] x0, double t1)
					throws InterruptedException {
		checkArgument(runs > 0, "Expected runs > 0");

		final List<TrajectoryRecorder> trList = Collections.synchronizedList(new LinkedList<TrajectoryRecorder>());

		final AtomicBoolean runFailed = new AtomicBoolean(false);
		final AtomicBoolean throwableSet = new AtomicBoolean(false);
		final Throwable[] throwable = new Throwable[1];
        final Semaphore jobQueueCounter = new Semaphore(maxQueuedJobs);
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
					if (runFailed.compareAndSet(false, true)) {
						throwable[0] = t;
						throwableSet.set(true);
					}
				}

				@Override
				public void onSuccess(TrajectoryRecorder tr) {
					jobQueueCounter.release();
//					releaseSubmission(queuedJobs);
					trList.add(tr);
				}

    		});

    		if (throwableSet.get())
    			throw new RuntimeException(throwable[0]);
        }

        executor.shutdown();
		while (!executor.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS));

        return trList;
	}

	@Override
	public FiniteStatisticalSummaryTrajectory simulateTrajectoryDistribution(
			int runs, ObjProvider<? extends T> modelProvider, ObjProvider<? extends FiniteTrajectoryRecorder> trProvider,
			double t0, double[] x0, double t1)
					throws InterruptedException {
		DummyTrajectoryMapper mapper = new DummyTrajectoryMapper();
		return simulateTrajectoryDistribution(runs, modelProvider, trProvider, mapper, t0, x0, t1);
	}

	@Override
	public FiniteStatisticalSummaryTrajectory simulateTrajectoryDistribution(
			int runs, ObjProvider<? extends T> modelProvider, ObjProvider<? extends FiniteTrajectoryRecorder> trProvider,
			final FiniteTrajectoryMapper mapper, double t0, double[] x0, double t1)
					throws InterruptedException {
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

        final Semaphore jobQueueCounter = new Semaphore(DEFAULT_MAX_QUEUED_JOBS);
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
		while (!executor.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS));
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
