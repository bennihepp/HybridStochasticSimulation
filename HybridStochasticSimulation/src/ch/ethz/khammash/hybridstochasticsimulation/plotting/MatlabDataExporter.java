package ch.ethz.khammash.hybridstochasticsimulation.plotting;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JOptionPane;

import com.jmatio.io.MatFileWriter;
import com.jmatio.types.MLArray;
import com.jmatio.types.MLChar;
import com.jmatio.types.MLDouble;
import com.jmatio.types.MLStructure;

public class MatlabDataExporter {

	public List<MLArray> buildMatlabData(List<PlotData> plotDataList) {
		return buildMatlabData(plotDataList, 2, Math.round(plotDataList.size() / 2));
	}
	public List<MLArray> buildMatlabData(List<PlotData> plotDataList, int rows, int cols) {
		int[] plotDim = { plotDataList.size(), 1 };
		MLStructure mlPlots = new MLStructure("plots", plotDim);
		for (int i=0; i < plotDataList.size(); i++ ) {
			PlotData plotData = plotDataList.get(i);
			int[] trajectoriesDim = { plotData.getNumberOfStates(), 1 };
	        MLStructure mlTrajectories = new MLStructure("trajectories", trajectoriesDim);
			if (plotData instanceof TrajectoryPlotData) {
				TrajectoryPlotData td = (TrajectoryPlotData)plotData;
				double[] tSeries = td.gettVector().toArray();
				MLDouble mltSeries = new MLDouble("tSeries", tSeries, tSeries.length);
				mlPlots.setField("tSeries", mltSeries, i);
				for (int j=0; j < plotData.getNumberOfStates(); j++) {
					double[] xSeries = td.getxVector(j).toArray();
					MLDouble mlxSeries = new MLDouble("xSeries", xSeries, xSeries.length);
					mlTrajectories.setField("xSeries", mlxSeries, j);
					String name = td.getName(j);
					mlTrajectories.setField("name", new MLChar("name", name), j);
				}
			} else if (plotData instanceof TrajectoryDistributionPlotData) {
				TrajectoryDistributionPlotData tdd = (TrajectoryDistributionPlotData)plotData;
				double[] tSeries = tdd.gettVector().toArray();
				MLDouble mltSeries = new MLDouble("tSeries", tSeries, tSeries.length);
				mlPlots.setField("tSeries", mltSeries, i);
				for (int j=0; j < plotData.getNumberOfStates(); j++) {
					double[] xSeries = tdd.getxVector(j).toArray();
					MLDouble mlxSeries = new MLDouble("xSeries", xSeries, xSeries.length);
					mlTrajectories.setField("xSeries", mlxSeries, j);
					double[] xSeriesStdDev = tdd.getxStdDevVector(j).toArray();
					MLDouble mlxSeriesStdDev = new MLDouble("xSeriesStdDev", xSeriesStdDev, xSeriesStdDev.length);
					mlTrajectories.setField("xSeriesStdDev", mlxSeriesStdDev, j);
					String name = tdd.getName(j);
					mlTrajectories.setField("name", new MLChar("name", name), j);
				}
			} else
				throw new UnsupportedOperationException("Unsupported implementation of PlotData");
			mlPlots.setField("trajectories", mlTrajectories, i);
			mlPlots.setField("title", new MLChar("title", plotData.getTitle()), i);
		}
		ArrayList<MLArray> list = new ArrayList<MLArray>(1);
		list.add(mlPlots);
		double[] rowArr = { rows };
		double[] colArr = { cols };
		list.add(new MLDouble("rows", rowArr, 1));
		list.add(new MLDouble("cols", colArr, 1));
		return list;
	}

	public void writeMatlabDataToFile(File file, List<MLArray> matlabData) {
		try {
			new MatFileWriter(file, matlabData);
		} catch (IOException e1) {
			String errorMsg = "Failed to export Matlab file:\n" + e1.toString();
			JOptionPane.showMessageDialog(null, errorMsg, "Error", JOptionPane.ERROR);
		}
	}

}
