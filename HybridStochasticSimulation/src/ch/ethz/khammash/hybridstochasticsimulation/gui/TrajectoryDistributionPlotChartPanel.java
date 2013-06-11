package ch.ethz.khammash.hybridstochasticsimulation.gui;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Paint;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.stat.descriptive.StatisticalSummary;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.LegendItem;
import org.jfree.chart.LegendItemCollection;
import org.jfree.chart.LegendItemSource;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import ch.ethz.khammash.hybridstochasticsimulation.matlab.MatlabPlotter;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.TrajectoryDistributionPlotData;

import com.google.common.collect.Iterators;

public class TrajectoryDistributionPlotChartPanel extends ChartPanel {

	private static final long serialVersionUID = -8825961199562105387L;

	final public Color[] DEFAULT_COLOR_ORDER;
	{
		DEFAULT_COLOR_ORDER = new Color[MatlabPlotter.DEFAULT_COLOR_ORDER.length];
		for (int i=0; i < DEFAULT_COLOR_ORDER.length; i++) {
			double[] color = MatlabPlotter.DEFAULT_COLOR_ORDER[i];
			DEFAULT_COLOR_ORDER[i] = new Color((float)color[0], (float)color[1], (float)color[2]);
		}
	}
//	final public Color[] DEFAULT_COLOR_ORDER = {
//			Color.blue,
//			Color.red,
//			Color.green,
//			Color.cyan,
//			Color.magenta,
//			Color.orange,
//			Color.pink,
//			Color.black,
//			Color.darkGray,
//			Color.lightGray,
//	};

	private final XYLineAndShapeRenderer renderer;
	private final XYSeriesCollection seriesCollection;
	private final BasicStroke stroke;
	private final BasicStroke dottedStroke;
	private List<String> speciesNames;
	private Color[] colorList;

	public TrajectoryDistributionPlotChartPanel() {
		super(ChartFactory.createXYLineChart(
				"Trajectory Distribution", // chart title
				"t", // x axis label
				"x", // y axis label
				new XYSeriesCollection(), // data
				PlotOrientation.VERTICAL, true, // include legend
				true, // tooltips
				false // urls
				));
		colorList = DEFAULT_COLOR_ORDER;
		seriesCollection = (XYSeriesCollection) getChart().getXYPlot().getDataset();
		JFreeChart chart = getChart();
		XYPlot plot = chart.getXYPlot();
		plot.setBackgroundPaint(Color.white);
		renderer = new XYLineAndShapeRenderer(true, false);
		plot.setRenderer(renderer);
		stroke = new BasicStroke(1.0f);
		dottedStroke = new BasicStroke(0.25f);
		speciesNames = new ArrayList<String>();
		LegendItemSource[] sources = {new LegendItemSource() {
			public LegendItemCollection getLegendItems() {
				LegendItemCollection lic = new LegendItemCollection();
				Iterator<Color> it = Iterators.cycle(colorList);
				for (int i = 0; i < speciesNames.size(); i++) {
					LegendItem li = new LegendItem(speciesNames.get(i), it.next());
					li.setShape(new Rectangle(5, 5));
					lic.add(li);
				}
				return lic;
			}
		}};
		chart.getLegend().setSources(sources);
	}

	public void setTitle(String title) {
		getChart().setTitle(title);
	}

	public void addSpecies(String name, XYSeries meanSeries, XYSeries stdDevPlusSeries, XYSeries stdDevMinusSeries) {
		int i = seriesCollection.getSeriesCount();
		seriesCollection.addSeries(meanSeries);
		seriesCollection.addSeries(stdDevPlusSeries);
		seriesCollection.addSeries(stdDevMinusSeries);
		speciesNames.add(name);
		for (int j=i; j < i + 3; j++) {
			Paint color = colorList[(j / 3) % colorList.length];
			if (j % 3 == 0)
				renderer.setSeriesStroke(j, stroke);
			else
				renderer.setSeriesStroke(j, dottedStroke);
			renderer.setSeriesPaint(j, color);
		}
	}

	public void addSpecies(String name, double[] tSeries, StatisticalSummary[] statistics, double plotScale) {
		XYSeries meanSeries = new XYSeries(name);
		XYSeries stdDevPlusSeries = new XYSeries(name + "+sigma");
		XYSeries stdDevMinusSeries = new XYSeries(name + "-sigma");
		for (int i = 0; i < statistics.length; i++) {
			double xMean = statistics[i].getMean();
			double xStdDev = statistics[i].getStandardDeviation();
			meanSeries.add(tSeries[i], plotScale * xMean);
			stdDevPlusSeries.add(tSeries[i], plotScale * (xMean + xStdDev));
			stdDevMinusSeries.add(tSeries[i], plotScale * (xMean - xStdDev));
		}
		addSpecies(name, meanSeries, stdDevPlusSeries, stdDevMinusSeries);
	}

	public void addSpecies(String[] names, double[] tSeries, StatisticalSummary[][] statistics, double[] plotScale) {
		for (int s=0; s < names.length; ++s) {
			StatisticalSummary[] tempStatistics = new StatisticalSummary[statistics.length];
			for (int i=0; i < statistics.length; i++)
				tempStatistics[i] = statistics[i][s];
			addSpecies(names[s], tSeries, tempStatistics, plotScale[s]);
		}
	}

	public void addSpecies(String name, RealVector tVector, RealVector meanVector, RealVector stdDevVector, double plotScale) {
		XYSeries meanSeries = new XYSeries(name);
		XYSeries stdDevPlusSeries = new XYSeries(name + "+sigma");
		XYSeries stdDevMinusSeries = new XYSeries(name + "-sigma");
		for (int i = 0; i < tVector.getDimension(); i++) {
			double time = tVector.getEntry(i);
			double xMean = meanVector.getEntry(i);
			double xStdDev = stdDevVector.getEntry(i);
			meanSeries.add(time, plotScale * xMean);
			stdDevPlusSeries.add(time, plotScale * (xMean + xStdDev));
			stdDevMinusSeries.add(time, plotScale * (xMean - xStdDev));
		}
		addSpecies(name, meanSeries, stdDevPlusSeries, stdDevMinusSeries);
	}

	public void addDistributionPlotData(TrajectoryDistributionPlotData tdd) {
		for (int s=0; s < tdd.getNumberOfStates(); ++s) {
			RealVector tVector = tdd.gettVector();
			RealVector xMeanVector = tdd.getxMeanVector(s);
			RealVector xStdDevVector = tdd.getxStdDevVector(s);
			addSpecies(tdd.getStateName(s), tVector, xMeanVector, xStdDevVector, tdd.getPlotScale(s));
		}
		setTitle(tdd.getTitle());
	}

}
