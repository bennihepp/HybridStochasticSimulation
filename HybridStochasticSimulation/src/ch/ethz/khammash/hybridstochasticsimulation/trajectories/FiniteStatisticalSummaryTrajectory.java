package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import org.apache.commons.math3.stat.descriptive.StatisticalSummary;

public interface FiniteStatisticalSummaryTrajectory extends FiniteDistributionTrajectory {

	StatisticalSummary[] getxStatisticalSummaryState(int timePoint);

	StatisticalSummary[] getxStatisticalSummarySeries(int state);

	StatisticalSummary[][] getxStatisticalSummary();

}
