package ch.ethz.bhepp.hybridstochasticsimulation.matlab;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;

import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.FinitePlotData;
import ch.ethz.bhepp.hybridstochasticsimulation.trajectories.VectorFiniteDistributionPlotData;

import com.jmatio.io.MatFileWriter;
import com.jmatio.types.MLArray;
import com.jmatio.types.MLCell;
import com.jmatio.types.MLChar;
import com.jmatio.types.MLDouble;
import com.jmatio.types.MLStructure;

public class MatlabDataExporter {

	public static MLArray buildString(String name, String value) {
		return new MLChar(name, value);
	}

	public static MLArray buildDouble(String name, double value) {
		double[] values = new double[1];
		return buildDouble(name, values);
	}

	public static MLDouble buildDouble(String name, double[] values) {
		return new MLDouble(name, values, values.length);
	}

	public static MLDouble buildDouble(String name, double[][] data) {
		return new MLDouble(name, data);
	}

	public static List<MLArray> buildMatlabPlotList(List<FinitePlotData> plotDataList) {
		MLStructure mlPlots = buildMatlabPlotData(plotDataList);
		List<MLArray> list = new LinkedList<MLArray>();
		list.add(mlPlots);
		return list;
	}

	public static MLStructure buildMatlabPlotData(List<FinitePlotData> plotDataList) {
		int[] plotDim = { plotDataList.size(), 1 };
		MLStructure mlPlots = new MLStructure("plots", plotDim);
		for (int i=0; i < plotDataList.size(); i++ ) {
			FinitePlotData plotData = plotDataList.get(i);
			buildMatlabPlotData(mlPlots, i, plotData);
		}
		return mlPlots;
	}

	public static MLArray buildMatlabPlotData(MLStructure mlPlots, int index, FinitePlotData plotData) {
		int[] trajectoriesDim = { plotData.getNumberOfStates(), 1 };
        MLStructure mlTrajectories = new MLStructure("trajectories", trajectoriesDim);
		double[] tSeries = plotData.gettVector().toArray();
		MLDouble mltSeries = buildDouble("tSeries", tSeries);
		mlPlots.setField("tSeries", mltSeries, index);
		Array2DRowRealMatrix xMatrix = new Array2DRowRealMatrix(tSeries.length, plotData.getNumberOfStates());
		MLDouble mlxMatrix = buildDouble("xMatrix", xMatrix.getData());
		mlPlots.setField("xMatrix", mlxMatrix, index);
		for (int j=0; j < plotData.getNumberOfStates(); j++) {
			double[] xSeries = plotData.getxVector(j).toArray();
			MLDouble mlxSeries = buildDouble("xSeries", xSeries);
			mlTrajectories.setField("xSeries", mlxSeries, j);
			if (plotData instanceof VectorFiniteDistributionPlotData) {
				VectorFiniteDistributionPlotData tdd = (VectorFiniteDistributionPlotData)plotData;
				double[] xSeriesStdDev = tdd.getxStdDevVector(j).toArray();
				MLDouble mlxSeriesStdDev = buildDouble("xSeriesStdDev", xSeriesStdDev);
				mlTrajectories.setField("xSeriesStdDev", mlxSeriesStdDev, j);
			}
			String name = plotData.getStateName(j);
			mlTrajectories.setField("name", buildString("name", name), j);
		}
		mlPlots.setField("trajectories", mlTrajectories, index);
		mlPlots.setField("title", buildString("title", plotData.getDescription()), index);
		return mlPlots;
	}

	public static List<MLArray> buildMatlabSimulationList(Map<String, List<FinitePlotData>> plotDataListMap) {
		int[] simulationDim = { plotDataListMap.size(), 1 };
//		MLStructure mlSimulations = new MLStructure("simulations", simulationDim);
		MLCell mlSimulations = new MLCell("simulations", simulationDim);
		int i = 0;
		for (Entry<String, List<FinitePlotData>> entry : plotDataListMap.entrySet()) {
			MLStructure mlPlots = buildMatlabPlotData(entry.getValue());
//			mlSimulations.setField("name", buildString("name", entry.getKey()), i);
//			mlSimulations.setField("plots", mlPlots, i);
			mlPlots.name = entry.getKey();
			mlSimulations.set(mlPlots, i);
			i++;
		}
		List<MLArray> list = new LinkedList<MLArray>();
		list.add(mlSimulations);
		return list;
	}

	public static void writeMatlabDataToFile(File file, Collection<MLArray> matlabData) throws IOException {
		new MatFileWriter(file, matlabData);
	}

	public static void writeMatlabDataToFile(File file, MLArray matlabData) throws IOException {
		MatFileWriter writer = new MatFileWriter();
		MLArray[] matlabDataArray = { matlabData };
		writer.write(file, Arrays.asList(matlabDataArray));
	}

}
