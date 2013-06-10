package ch.ethz.khammash.hybridstochasticsimulation.plotting;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.util.FastMath;

import com.jmatio.io.MatFileWriter;
import com.jmatio.types.MLArray;
import com.jmatio.types.MLChar;
import com.jmatio.types.MLDouble;
import com.jmatio.types.MLStructure;

public class MatlabDataExporter {

	public List<MLArray> buildMatlabData(List<PlotData> plotDataList) {
		return buildMatlabData(plotDataList, 2, FastMath.round(plotDataList.size() / 2));
	}

	public List<MLArray> buildMatlabData(List<PlotData> plotDataList, int rows, int cols) {
		int[] plotDim = { plotDataList.size(), 1 };
		MLStructure mlPlots = new MLStructure("plots", plotDim);
		for (int i=0; i < plotDataList.size(); i++ ) {
			PlotData plotData = plotDataList.get(i);
			int[] trajectoriesDim = { plotData.getNumberOfStates(), 1 };
	        MLStructure mlTrajectories = new MLStructure("trajectories", trajectoriesDim);
			double[] tSeries = plotData.gettVector().toArray();
			MLDouble mltSeries = new MLDouble("tSeries", tSeries, tSeries.length);
			mlPlots.setField("tSeries", mltSeries, i);
			Array2DRowRealMatrix xMatrix = new Array2DRowRealMatrix(tSeries.length, plotData.getNumberOfStates());
			MLDouble mlxMatrix = new MLDouble("xMatrix", xMatrix.getData());
			mlPlots.setField("xMatrix", mlxMatrix, i);
			for (int j=0; j < plotData.getNumberOfStates(); j++) {
				double[] xSeries = plotData.getxVector(j).toArray();
				MLDouble mlxSeries = new MLDouble("xSeries", xSeries, xSeries.length);
				mlTrajectories.setField("xSeries", mlxSeries, j);
				if (plotData instanceof TrajectoryDistributionPlotData) {
					TrajectoryDistributionPlotData tdd = (TrajectoryDistributionPlotData)plotData;
					double[] xSeriesStdDev = tdd.getxStdDevVector(j).toArray();
					MLDouble mlxSeriesStdDev = new MLDouble("xSeriesStdDev", xSeriesStdDev, xSeriesStdDev.length);
					mlTrajectories.setField("xSeriesStdDev", mlxSeriesStdDev, j);
				}
				String name = plotData.getStateName(j);
				mlTrajectories.setField("name", new MLChar("name", name), j);
			}
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

	public void writeMatlabDataToFile(File file, List<MLArray> matlabData) throws IOException {
		new MatFileWriter(file, matlabData);
	}

}
