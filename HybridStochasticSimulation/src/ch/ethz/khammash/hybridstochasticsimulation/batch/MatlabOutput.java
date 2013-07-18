package ch.ethz.khammash.hybridstochasticsimulation.batch;

import static com.google.common.base.Preconditions.checkArgument;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import ch.ethz.khammash.hybridstochasticsimulation.matlab.MatlabDataExporter;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FinitePlotData;

import com.jmatio.types.MLArray;


public class MatlabOutput implements SimulationOutput {

    private Map<String, List<FinitePlotData>> plotDataListMap;
    private File outputFile;
    private int rows = -1;
    private int cols = -1;

    public MatlabOutput(String outputFilename, boolean overwrite) throws IOException {
        this(new File(outputFilename), overwrite);
    }

    public MatlabOutput(File outputFile, boolean overwrite) throws IOException {
        plotDataListMap = new HashMap<>();
        this.outputFile = outputFile;
        if (!overwrite)
            checkArgument(!outputFile.exists(), "Output file already exists!");
        try {
            outputFile.createNewFile();
        } catch (IOException e) {
            System.err.println("Couldn't create file for output " + outputFile.getAbsolutePath());
            throw e;
        }
        checkArgument(outputFile.canWrite());
    }

    @Override
    public String toString() {
        return outputFile.getAbsolutePath();
    }

    @Override
    public void add(String simulationName, FinitePlotData plotData) {
        getPlotDataList(simulationName).add(plotData);
    }

    @Override
    public void addAll(String simulationName, List<FinitePlotData> plotDataList) {
        getPlotDataList(simulationName).addAll(plotDataList);
    }

    private List<FinitePlotData> getPlotDataList(String simulationName) {
        if (!plotDataListMap.containsKey(simulationName))
            plotDataListMap.put(simulationName, new LinkedList<FinitePlotData>());
        return plotDataListMap.get(simulationName);
    }

    @Override
    public void write() throws IOException {
        MatlabDataExporter mde = new MatlabDataExporter();
        List<MLArray> matlabData = mde.buildMatlabSimulationList(plotDataListMap);
        if (getRows() > 0)
            matlabData.add(mde.buildDouble("rows", getRows()));
        if (getCols() > 0)
            matlabData.add(mde.buildDouble("cols", getCols()));
        mde.writeMatlabDataToFile(outputFile, matlabData);
    }

    public int getRows() {
        return rows;
    }

    public void setRows(int rows) {
        this.rows = rows;
    }

    public int getCols() {
        return cols;
    }

    public void setCols(int cols) {
        this.cols = cols;
    }

}
