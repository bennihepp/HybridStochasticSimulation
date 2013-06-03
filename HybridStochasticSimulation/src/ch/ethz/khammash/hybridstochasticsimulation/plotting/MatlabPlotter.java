package ch.ethz.khammash.hybridstochasticsimulation.plotting;

import java.util.List;

import matlabcontrol.MatlabInvocationException;
import matlabcontrol.MatlabProxy;
import matlabcontrol.extensions.MatlabNumericArray;
import matlabcontrol.extensions.MatlabTypeConverter;

public class MatlabPlotter {

	final double[][] colorOrder = {
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

	public void plot(List<PlotData> plotDataList, int rows, int cols) throws MatlabInvocationException {
		MatlabProxy proxy = session.getProxy();

		if (session.isExistingSession()) {
			proxy.eval("clear java");
			proxy.eval("clear all");
		}

		MatlabTypeConverter processor = new MatlabTypeConverter(proxy);
	    processor.setNumericArray("colorOrder", new MatlabNumericArray(colorOrder, null));

	    proxy.eval("figure");
	    proxy.eval("clf");
	    proxy.eval("orient landscape");
	    proxy.eval("set(gcf, 'PaperUnits', 'centimeters')");
	    proxy.eval("set(gcf, 'PaperType', 'A4')");
	    for (int i=0; i < plotDataList.size(); i++) {
	    	PlotData plotData = plotDataList.get(i);
		    proxy.eval("subplot(" + rows + "," + cols + "," + (i+1) + ")");
		    proxy.eval("hold on");
			String[] plotLabels = new String[plotData.getNumberOfStates()];
			if (plotData instanceof TrajectoryPlotData) {
				TrajectoryPlotData td = (TrajectoryPlotData)plotData;
				double[] tSeries = td.gettVector().toArray();
				proxy.setVariable("tSeries", tSeries);
				for (int j=0; j < plotData.getNumberOfStates(); j++) {
					double[] xSeries = td.getxVector(j).toArray();
					proxy.setVariable("xSeries", xSeries);
					int colorOrderIndex = (j % colorOrder.length) + 1;
					proxy.eval("plot(tSeries, xSeries, 'Color', colorOrder("+colorOrderIndex+",:))");
					plotLabels[j] = td.getName(j);
				}
			} else if (plotData instanceof TrajectoryDistributionPlotData) {
				TrajectoryDistributionPlotData tdd = (TrajectoryDistributionPlotData)plotData;
				double[] tSeries = tdd.gettVector().toArray();
				proxy.setVariable("tSeries", tSeries);
				for (int j=0; j < plotData.getNumberOfStates(); j++) {
					double[] xSeries = tdd.getxVector(j).toArray();
					proxy.setVariable("xSeries", xSeries);
					int colorOrderIndex = (j % colorOrder.length) + 1;
					proxy.eval("plot(tSeries, xSeries, 'Color', colorOrder("+colorOrderIndex+",:))");
					double[] xSeriesStdDev = tdd.getxStdDevVector(j).toArray();
					proxy.setVariable("xSeriesStdDev", xSeriesStdDev);
					proxy.eval("plot(tSeries, xSeries+xSeriesStdDev, 'Color', colorOrder("+colorOrderIndex+",:))");
					proxy.eval("plot(tSeries, xSeries-xSeriesStdDev, 'Color', colorOrder("+colorOrderIndex+",:))");
					plotLabels[j] = tdd.getName(j);
				}
			} else
				throw new UnsupportedOperationException("Unsupported implementation of PlotData");
			proxy.setVariable("plotLabels", plotLabels);
			proxy.eval("legend(plotLabels)");
			proxy.setVariable("plotTitle", plotData.getTitle());
			proxy.eval("title(plotTitle)");
			proxy.setVariable("xAxisLabel", "Time in seconds");
			proxy.setVariable("yAxisLabel", "Copy number");
			proxy.eval("xlabel(xAxisLabel)");
			proxy.eval("ylabel(yAxisLabel)");
		    proxy.eval("hold off");
	    }
	}

}
