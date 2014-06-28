package ch.ethz.bhepp.hybridstochasticsimulation.matlab;

import java.util.ArrayList;
import java.util.List;

import matlabcontrol.MatlabInvocationException;
import matlabcontrol.MatlabProxy;
import matlabcontrol.extensions.MatlabNumericArray;
import matlabcontrol.extensions.MatlabTypeConverter;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FinitePlotData;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.VectorFiniteDistributionPlotData;

public class MatlabPlotter {

	final public static double[][] DEFAULT_COLOR_ORDER = {
		      {     0,         0,    1.0000},
		      {1.0000,         0,         0},
		      {     0,    1.0000,         0},
		      {     0,         0,    0.1724},
		      {1.0000,    0.1034,    0.7241},
		      {1.0000,    0.8276,         0},
		      {     0,    0.3448,         0},
		      {0.5172,    0.5172,    1.0000},
		      {0.6207,    0.3103,    0.2759},
		      {     0,    1.0000,    0.7586},
		      {     0,    0.5172,    0.5862},
		      {     0,         0,    0.4828},
		      {0.5862,    0.8276,    0.3103},
		      {0.9655,    0.6207,    0.8621},
		      {0.8276,    0.0690,    1.0000},
		      {0.4828,    0.1034,    0.4138},
		      {0.9655,    0.0690,    0.3793},
		      {1.0000,    0.7586,    0.5172},
		      {0.1379,    0.1379,    0.0345},
		      {0.5517,    0.6552,    0.4828},
	};

	MatlabSession session;

	public MatlabPlotter(MatlabSession session) {
		this.session = session;
	}

	public void plot(FinitePlotData plotData) throws MatlabInvocationException {
		List<FinitePlotData> plotDataList = new ArrayList<>(1);
		plotDataList.add(plotData);
		plot(plotDataList, 1, 1);
	}

	public void plot(List<FinitePlotData> plotDataList, int rows, int cols) throws MatlabInvocationException {
		MatlabProxy proxy = session.getProxy();

		if (session.isExistingSession()) {
			proxy.eval("clear java");
			proxy.eval("clear all");
		}

		boolean distinguishableColorFunctionExists = false;
		double[] result = (double[])proxy.returningEval("exist('distinguishable_colors')", 1)[0];
		if (result[0] == 2.0)
			distinguishableColorFunctionExists = true;
		else {
			MatlabTypeConverter processor = new MatlabTypeConverter(proxy);
		    processor.setNumericArray("colorOrder", new MatlabNumericArray(DEFAULT_COLOR_ORDER, null));
		}

	    proxy.eval("figure");
	    proxy.eval("clf");
	    proxy.eval("orient landscape");
	    proxy.eval("set(gcf, 'PaperUnits', 'centimeters')");
	    proxy.eval("set(gcf, 'PaperType', 'A4')");
	    for (int i=0; i < plotDataList.size(); i++) {
	    	FinitePlotData plotData = plotDataList.get(i);
		    proxy.eval("subplot(" + rows + "," + cols + "," + (i+1) + ")");
		    proxy.eval("hold on");
		    if (distinguishableColorFunctionExists)
				proxy.eval("colorOrder = distinguishable_colors(" + plotData.getNumberOfStates() + ")");
			String[] plotLabels = new String[plotData.getNumberOfStates()];
			double[] tSeries = plotData.gettVector().toArray();
			proxy.setVariable("tSeries", tSeries);
			proxy.eval("plotHandles = zeros(1," + plotData.getNumberOfStates() + ")");
			for (int j=0; j < plotData.getNumberOfStates(); j++) {
				double[] xSeries = plotData.getxVector(j).toArray();
				proxy.setVariable("xSeries", xSeries);
				String colorOrderExpression;
				if (!distinguishableColorFunctionExists) {
					int colorOrderIndex = (j % DEFAULT_COLOR_ORDER.length) + 1;
					colorOrderExpression = "colorOrder(" + colorOrderIndex + ",:)";
				} else
					colorOrderExpression = "colorOrder(" + (j+1) + ",:)";
				proxy.eval("plotHandles(" + (j+1) + ") = plot(tSeries, xSeries, 'Color', " + colorOrderExpression + ")");
				if (plotData instanceof VectorFiniteDistributionPlotData) {
					VectorFiniteDistributionPlotData tdd = (VectorFiniteDistributionPlotData)plotData;
					double[] xSeriesStdDev = tdd.getxStdDevVector(j).toArray();
					proxy.setVariable("xSeriesStdDev", xSeriesStdDev);
					proxy.eval("plot(tSeries, xSeries+xSeriesStdDev, ':', 'Color', " + colorOrderExpression + ")");
					proxy.eval("plot(tSeries, xSeries-xSeriesStdDev, ':', 'Color', " + colorOrderExpression + ")");
				}
				plotLabels[j] = plotData.getStateName(j);
			}
			proxy.setVariable("plotLabels", plotLabels);
			proxy.eval("legend(plotHandles, plotLabels, 'Location', 'Best')");
			proxy.setVariable("plotTitle", plotData.getDescription());
			proxy.eval("title(plotTitle)");
			proxy.setVariable("xAxisLabel", "Time in seconds");
			proxy.setVariable("yAxisLabel", "Copy number");
			proxy.eval("xlabel(xAxisLabel)");
			proxy.eval("ylabel(yAxisLabel)");
		    proxy.eval("hold off");
	    }
	}

}
