package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.concurrent.atomic.AtomicInteger;

import org.apache.commons.math3.stat.descriptive.SynchronizedSummaryStatistics;

public class VectorFiniteDistributionTrajectoryBuilder implements FiniteDistributionTrajectoryBuilder {

	private volatile boolean initialized = false;
	private double[] tSeries;
	private SynchronizedSummaryStatistics[][] xSeriesStatistics;
	private AtomicInteger numberOfAddedTrajectories = new AtomicInteger(0);

	@Override
	public void addTrajectory(FiniteTrajectory tr) {
		// Double-locking idiom (this is important for thread-safety!
		// Don't change this unless you really know what you are doing!
		if (initialized == false) {
			synchronized (this) {
				if (initialized == false) {
					tSeries = tr.gettSeries();
					xSeriesStatistics = new SynchronizedSummaryStatistics[tr.getNumberOfStates()][tr.getNumberOfTimePoints()];
			        for (int s=0; s < tr.getNumberOfStates(); s++)
			        	for (int i=0; i < tr.getNumberOfTimePoints(); i++)
							xSeriesStatistics[s][i] = new SynchronizedSummaryStatistics();
					initialized = true;
				}
			}
		}
		// Check that the dimensions of the trajectory are consistent.
		checkArgument(tr.getNumberOfStates() == xSeriesStatistics.length);
		if (xSeriesStatistics.length > 0)
			checkArgument(tr.getNumberOfTimePoints() == xSeriesStatistics[0].length);
		numberOfAddedTrajectories.incrementAndGet();
		// Add knew values to the statistical summary objects
		for (int s = 0; s < tr.getNumberOfStates(); s++) {
			for (int i = 0; i < tr.getNumberOfTimePoints(); i++) {
				xSeriesStatistics[s][i].addValue(tr.getxState(s, i));
			}
		}
	}

	@Override
	public FiniteStatisticalSummaryTrajectory getDistributionTrajectory() {
		return VectorFiniteStatisticalSummaryPlotData.createFromStatisticalSummary(tSeries, xSeriesStatistics);
	}

	@Override
	public int getNumberOfAddedTrajectories() {
		return numberOfAddedTrajectories.get();
	}

}
