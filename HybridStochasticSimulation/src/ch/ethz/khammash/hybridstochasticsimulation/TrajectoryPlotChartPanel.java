package ch.ethz.khammash.hybridstochasticsimulation;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Paint;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.LegendItem;
import org.jfree.chart.LegendItemCollection;
import org.jfree.chart.LegendItemSource;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.renderer.xy.XYStepRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class TrajectoryPlotChartPanel extends ChartPanel {

	private static final long serialVersionUID = -8825961199562105387L;

	XYLineAndShapeRenderer renderer;
	final XYSeriesCollection seriesCollection;
	final BasicStroke stroke;
	final Color[] colorList = { Color.blue, Color.red, Color.green, Color.cyan, Color.magenta, Color.orange };
	List<String> speciesNames;

	public TrajectoryPlotChartPanel() {
		this(false);
	}

	public TrajectoryPlotChartPanel(boolean useStepRenderer) {
		super(ChartFactory.createXYLineChart(
				"Trajectory Distribution", // chart title
				"t", // x axis label
				"x", // y axis label
				new XYSeriesCollection(), // data
				PlotOrientation.VERTICAL, true, // include legend
				true, // tooltips
				false // urls
				));
		seriesCollection = (XYSeriesCollection) getChart().getXYPlot().getDataset();
		JFreeChart chart = getChart();
		XYPlot plot = chart.getXYPlot();
		plot.setBackgroundPaint(Color.white);
		if (useStepRenderer)
			renderer = new XYStepRenderer();
		else
			renderer = new XYLineAndShapeRenderer(true, false);
		plot.setRenderer(renderer);
		stroke = new BasicStroke(1.0f);
		speciesNames = new ArrayList<String>();
		LegendItemSource[] sources = {new LegendItemSource() {
			public LegendItemCollection getLegendItems() {
				LegendItemCollection lic = new LegendItemCollection();
				for (int i = 0; i < speciesNames.size(); i++) {
					LegendItem li = new LegendItem(speciesNames.get(i), colorList[i]);
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

	public void setStepRenderer() {
		renderer = new XYStepRenderer();
		XYPlot plot = getChart().getXYPlot();
		plot.setRenderer(renderer);
	}

	public void setLineRenderer() {
		renderer = new XYLineAndShapeRenderer(true, false);
		XYPlot plot = getChart().getXYPlot();
		plot.setRenderer(renderer);
	}

	public void addSpecies(String name, double[] tSeries, double[] xSeries, double plotScale) {
		XYSeries xySeries = new XYSeries(name);
		for (int i = 0; i < xSeries.length; i++)
			xySeries.add(tSeries[i], plotScale * xSeries[i]);
		int i = seriesCollection.getSeriesCount();
		seriesCollection.addSeries(xySeries);
		speciesNames.add(name);
		Paint color = colorList[i % colorList.length];
		renderer.setSeriesStroke(i, stroke);
		renderer.setSeriesPaint(i, color);
	}

	public void addSpecies(String[] names, double[] tSeries, double[][] xSeries, double[] plotScale) {
		for (int s=0; s < names.length; ++s) {
			double[] tempXSeries = new double[xSeries.length];
			for (int i=0; i < xSeries.length; i++)
				tempXSeries[i] = xSeries[i][s];
			addSpecies(names[s], tSeries, tempXSeries, plotScale[s]);
		}
	}

	public void addSpecies(String name, RealVector tVector, RealVector xVector, double plotScale) {
		addSpecies(name, tVector.toArray(), xVector.toArray(), plotScale);
	}

	public void addSpecies(String[] names, RealVector tVector, RealMatrix xMatrix, double[] plotScale) {
		for (int s=0; s < names.length; ++s) {
			RealVector xSeries = xMatrix.getRowVector(s);
			addSpecies(names[s], tVector, xSeries, plotScale[s]);
		}
	}

	public void addSpecies(String[] names, TrajectoryPlotData td, double[] plotScale) {
		for (int s=0; s < names.length; ++s) {
			addSpecies(names[s], td.gettVector(), td.getxVector(s), plotScale[s]);
		}
	}

}
