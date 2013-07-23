package ch.ethz.khammash.hybridstochasticsimulation.trajectories;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.stat.descriptive.StatisticalSummary;

public class VectorFiniteStatisticalSummaryPlotData extends VectorFiniteDistributionPlotData implements FiniteStatisticalSummaryTrajectory {

	private static final long serialVersionUID = 988659663487743261L;

	private StatisticalSummary[][] xSeriesStatistics;

	public static FiniteStatisticalSummaryTrajectory createFromStatisticalSummary(
			double[] tSeries, StatisticalSummary[][] xSeriesStatistics) {
		RealVector tVector = new ArrayRealVector(tSeries);
		RealMatrix xMeanMatrix = new Array2DRowRealMatrix(xSeriesStatistics.length, xSeriesStatistics[0].length);
		RealMatrix xStdDevMatrix = xMeanMatrix.copy();
		for (int s = 0; s < xSeriesStatistics.length; s++)
			for (int i = 0; i < xSeriesStatistics[0].length; i++) {
				double xMean = xSeriesStatistics[s][i].getMean();
				double xStdDev = xSeriesStatistics[s][i].getStandardDeviation();
				xMeanMatrix.setEntry(s, i, xMean);
				xStdDevMatrix.setEntry(s, i, xStdDev);
			}
		return new VectorFiniteStatisticalSummaryPlotData(tVector, xMeanMatrix, xStdDevMatrix);
	}

	protected VectorFiniteStatisticalSummaryPlotData(RealVector tVector, RealMatrix xMeanMatrix, RealMatrix xStdDevMatrix) {
		super(tVector, xMeanMatrix, xStdDevMatrix);
	}

	protected void setxStatisticalSummary(StatisticalSummary[][] xSeriesStatistics) {
		this.xSeriesStatistics = xSeriesStatistics;
	}

	@Override
	public StatisticalSummary[] getxStatisticalSummarySeries(int state) {
		return xSeriesStatistics[state];
	}

	@Override
	public StatisticalSummary[][] getxStatisticalSummary() {
		return xSeriesStatistics;
	}

	@Override
	public StatisticalSummary[] getxStatisticalSummaryState(int timePoint) {
		StatisticalSummary[] stateStatistics = new StatisticalSummary[getNumberOfStates()];
		for (int s=0; s < getNumberOfStates(); s++)
			stateStatistics[s] = xSeriesStatistics[s][timePoint];
		return stateStatistics;
	}

}
