package ch.ethz.bhepp.hybridstochasticsimulation.controllers;

import org.apache.commons.math3.stat.descriptive.SynchronizedSummaryStatistics;

import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FiniteTrajectoryRecorder;

public class DefaultDistributionSimulationWorker implements DistributionSimulationWorker {

	private FiniteSimulationWorker worker;
	private double[] tSeries;
	private SynchronizedSummaryStatistics[][] xSeriesStatistics;

	DefaultDistributionSimulationWorker(FiniteSimulationWorker worker,
			double[] tSeries, SynchronizedSummaryStatistics[][] xSeriesStatistics) {
		this.worker = worker;
		this.tSeries = tSeries;
		this.xSeriesStatistics = xSeriesStatistics;
	}

	@Override
	public void run() {
		FiniteTrajectoryRecorder tr = worker.simulate();
		double[] tSeries = tr.gettSeries();
		for (int i = 0; i < tSeries.length; i++) {
			this.tSeries[i] = tSeries[i];
			double[] x = tr.getInterpolatedState(tSeries[i]);
			for (int s = 0; s < tr.getNumberOfStates(); s++) {
				xSeriesStatistics[s][i].addValue(x[s]);
			}
		}
	}

}
