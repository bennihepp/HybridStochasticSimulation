package ch.ethz.bhepp.hybridstochasticsimulation.trajectories;

import org.apache.commons.math3.stat.descriptive.StatisticalSummary;

public interface FiniteStatisticalSummaryTrajectory extends FiniteDistributionTrajectory {

	StatisticalSummary[] getxStatisticalSummaryState(int timePoint);

	StatisticalSummary[] getxStatisticalSummarySeries(int state);

	StatisticalSummary[][] getxStatisticalSummary();

}
