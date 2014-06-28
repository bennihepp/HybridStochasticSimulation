package ch.ethz.bhepp.hybridstochasticsimulation.gui;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Paint;
import java.util.ArrayList;
import java.util.Iterator;
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
import org.jfree.chart.renderer.xy.AbstractXYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.renderer.xy.XYStepRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import ch.ethz.bhepp.hybridstochasticsimulation.matlab.MatlabPlotter;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FinitePlotData;

import com.google.common.collect.Iterators;

public class TrajectoryPlotChartPanel extends ChartPanel {

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

	private AbstractXYItemRenderer stepRenderer;
	private AbstractXYItemRenderer continuousRenderer;
	private final XYSeriesCollection seriesCollection;
	private final BasicStroke stroke;
	private List<String> speciesNames;
	private Color[] colorList;

	public TrajectoryPlotChartPanel() {
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
		stepRenderer = new XYStepRenderer();
		continuousRenderer = new XYLineAndShapeRenderer(true, false);
		stroke = new BasicStroke(1.0f);
		speciesNames = new ArrayList<String>();
		LegendItemSource[] sources = {new LegendItemSource() {
			public LegendItemCollection getLegendItems() {
				LegendItemCollection lic = new LegendItemCollection();
				Iterator<Color> it = Iterators.cycle(colorList);
				for (int i = 0; i < speciesNames.size(); i++) {
					LegendItem li = new LegendItem(speciesNames.get(i), it.next());
//					li.setShape(new Rectangle(5, 5));
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

	public void addSpecies(String name, double[] tSeries, double[] xSeries, double plotScale, boolean isDiscrete) {
		XYSeries xySeries = new XYSeries(name);
		for (int i = 0; i < xSeries.length; i++)
			xySeries.add(tSeries[i], plotScale * xSeries[i]);
		int i = seriesCollection.getSeriesCount();
		seriesCollection.addSeries(xySeries);
		speciesNames.add(name);
		Paint color = colorList[i % colorList.length]; 
		AbstractXYItemRenderer renderer;
		if (isDiscrete)
			renderer = stepRenderer;
		else
			renderer = continuousRenderer;
		getChart().getXYPlot().setRenderer(i, renderer);
		renderer.setSeriesStroke(i, stroke);
		renderer.setSeriesPaint(i, color);
	}

	public void addSpecies(String[] names, double[] tSeries, double[][] xSeries, double[] plotScale, boolean isDiscrete) {
		for (int s=0; s < names.length; ++s) {
			double[] tempXSeries = new double[xSeries.length];
			for (int i=0; i < xSeries.length; i++)
				tempXSeries[i] = xSeries[i][s];
			addSpecies(names[s], tSeries, tempXSeries, plotScale[s], isDiscrete);
		}
	}

	public void addSpecies(String name, RealVector tVector, RealVector xVector, double plotScale, boolean isDiscrete) {
		addSpecies(name, tVector.toArray(), xVector.toArray(), plotScale, isDiscrete);
	}

	public void addSpecies(String[] names, RealVector tVector, RealMatrix xMatrix, double[] plotScale, boolean isDiscrete) {
		for (int s=0; s < names.length; ++s) {
			RealVector xSeries = xMatrix.getRowVector(s);
			addSpecies(names[s], tVector, xSeries, plotScale[s], isDiscrete);
		}
	}

	public void addPlotData(FinitePlotData pd) {
		for (int s=0; s < pd.getNumberOfStates(); ++s) {
			addSpecies(pd.getStateName(s), pd.gettVector(), pd.getxVector(s), pd.getPlotScale(s), pd.isDiscrete());
		}
		setTitle(pd.getDescription());
	}

}
