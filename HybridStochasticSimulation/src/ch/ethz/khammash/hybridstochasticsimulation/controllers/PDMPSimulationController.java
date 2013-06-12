package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import java.util.concurrent.ExecutorService;

import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.ContinuousTrajectoryRecorder;


public class PDMPSimulationController<T extends PDMPModel>
		extends AbstractSimulationController<T, ContinuousTrajectoryRecorder<T>> {

    public PDMPSimulationController() {
    	super();
		construct();
    }

	public PDMPSimulationController(int numOfThreads) {
		super(numOfThreads);
		construct();
    }

	public PDMPSimulationController(ExecutorService executor) {
		super(executor);
		construct();
    }

	final private void construct() {
		useDefaultPDMPSimulatorFactory();
	}

	final public void useDefaultPDMPSimulatorFactory() {
		useDefaultPDMPSimulatorFactory(new DefaultIntegratorFactory());
	}

	final public void useDefaultPDMPSimulatorFactory(IntegratorFactory integratorFactory) {
		setSimulatorFactory(new DefaultPDMPSimulatorFactory<T>(integratorFactory));
	}

//	public ContinuousPDMPTrajectory<PDMPModel> simulateTrajectory(PDMPModel model, double t0, double[] x0, double t1) {
//		PDMPSimulator simulator = createSimulator();
//		SimulationWorker<PDMPModel, ContinuousTrajectoryRecorder<PDMPModel>> sw
//			= new SimulationWorker(simulator, model, t0, x0, t1);
//		return sw.simulate();
//	}

//	public double[][] simulateTrajectory(PDMPModel model, double[] tSeries, double[] x0) {
//		FinitePDMPTrajectory<PDMPModel> mt = new FinitePDMPTrajectory<PDMPModel>(tSeries);
//		simulateTrajectory(model, mt, tSeries[0], x0, tSeries[tSeries.length - 1]);
//		return mt.getxSeries();
//	}

//	public void simulateTrajectory(PDMPModel model, FinitePDMPTrajectory<PDMPModel> mt, double t0, double[] x0, double t1) {
//		PDMPSimulator simulator = createSimulator();
//		FiniteSimulationWorker sw = new FiniteSimulationWorker(simulator, mt, model, t0, x0, t1);
//		sw.simulate();
//	}

//	public StatisticalSummary[][] simulateTrajectoryDistribution(int runs,
//			PDMPModelFactory modelFactory, double[] tSeries, double[] x0)
//			throws InterruptedException, CancellationException, ExecutionException {
//		SynchronizedSummaryStatistics[][] xSeriesStatistics = new SynchronizedSummaryStatistics[tSeries.length][x0.length];
//		for (int i = 0; i < tSeries.length; i++)
//			for (int s = 0; s < x0.length; s++)
//				xSeriesStatistics[i][s] = new SynchronizedSummaryStatistics();
//		for (int k = 0; k < runs; k++) {
//			PDMPSimulator simulator = createSimulator();
//			// PDMPModel has an internal state so we create a copy for each run
//			PDMPModel model = modelFactory.createModel();
//			Runnable r = new TrajectoryDistributionSimulator(simulator, model, xSeriesStatistics, tSeries, x0);
//			executor.execute(r);
//		}
//		executor.shutdown();
//		try {
//			executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
//		} catch (InterruptedException e) {
//			executor.shutdownNow();
//			Thread.currentThread().interrupt();
//		}
//		return xSeriesStatistics;
//	}

//	public StatisticalSummary[][] simulateFiniteTrajectoryDistribution(int runs,
//			PDMPModelFactory modelFactory,
//			FinitePDMPModelTrajectoryFactory mtFactory, double t0, double[] x0,
//			double t1) throws InterruptedException, CancellationException,
//			ExecutionException {
//		SynchronizedSummaryStatistics[][] xSeriesStatistics = null;
//		// We also want to have different random number sequences for each run
//		for (int k = 0; k < runs; k++) {
//			PDMPSimulator simulator = createSimulator();
//			// PDMPModel has an internal state so we create a copy for each run
//			PDMPModel model = modelFactory.createModel();
//			FinitePDMPTrajectory mt = mtFactory.createModelTrajectory();
//			if (xSeriesStatistics == null) {
//				xSeriesStatistics = new SynchronizedSummaryStatistics[x0.length][mt.gettSeries().length];
//				for (int s = 0; s < xSeriesStatistics.length; s++)
//					for (int i = 0; i < xSeriesStatistics[s].length; i++)
//						xSeriesStatistics[s][i] = new SynchronizedSummaryStatistics();
//			}
//			Runnable r = new FiniteTrajectoryDistributionSimulator(simulator,
//					mt, model, xSeriesStatistics, t0, x0.clone(), t1);
//			executor.execute(r);
//		}
//		executor.shutdown();
//		try {
//			executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
//		} catch (InterruptedException e) {
//			executor.shutdownNow();
//			Thread.currentThread().interrupt();
//		}
//		return xSeriesStatistics;
//	}

}
