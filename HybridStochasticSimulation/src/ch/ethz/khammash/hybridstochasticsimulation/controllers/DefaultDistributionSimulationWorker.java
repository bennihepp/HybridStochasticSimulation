package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import static com.google.common.base.Preconditions.checkArgument;

import org.apache.commons.math3.stat.descriptive.SynchronizedSummaryStatistics;

import ch.ethz.khammash.hybridstochasticsimulation.models.ReactionNetworkModel;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryRecorder;

public class DefaultDistributionSimulationWorker<T extends ReactionNetworkModel, E extends TrajectoryRecorder<T>>
		implements DistributionSimulationWorker {

	private SimulationWorker<T, E> worker;
	private double[] tSeries;
	private SynchronizedSummaryStatistics[][] xSeriesStatistics;

	DefaultDistributionSimulationWorker(SimulationWorker<T, E> worker,
			SynchronizedSummaryStatistics[][] xSeriesStatistics, double[] tSeries, double[] x0) {
		checkArgument(tSeries.length >= 2, "Expected tSeries.length >= 2");
		checkArgument(x0.length == xSeriesStatistics.length, "Expected x0.length == xStatistics.length");
		checkArgument(tSeries.length == xSeriesStatistics[0].length, "Expected tSeries.length == xStatistics[0].length");
		this.worker = worker;
		this.xSeriesStatistics = xSeriesStatistics;
		this.tSeries = tSeries;
	}

	@Override
	public void run() {
		E tr = worker.simulate();
		for (int i = 0; i < tSeries.length; i++) {
			double[] x = tr.getInterpolatedState(tSeries[i]);
			for (int s = 0; s < x.length; s++) {
				xSeriesStatistics[s][i].addValue(x[s]);
			}
		}
	}

}
