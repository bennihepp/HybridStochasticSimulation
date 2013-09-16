package ch.ethz.khammash.hybridstochasticsimulation.io;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob.OutputAlreadyWrittenException;
import ch.ethz.khammash.hybridstochasticsimulation.matlab.MatlabDataExporter;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FinitePlotData;

import com.jmatio.types.MLArray;


public class MatlabOutput implements SimulationOutput {

    private Map<String, List<FinitePlotData>> plotDataListMap;
    private File outputFile;
//    private int rows = -1;
//    private int cols = -1;
    private boolean outputWritten = false;
	private boolean overwrite;

    public MatlabOutput(String outputFilename, boolean overwrite) throws IOException {
        this(new File(outputFilename), overwrite);
    }

    public MatlabOutput(File outputFile, boolean overwrite) throws IOException {
        plotDataListMap = new HashMap<>();
        this.outputFile = outputFile;
        this.overwrite = overwrite;
    }

    @Override
    public String toString() {
        return outputFile.getAbsolutePath();
    }

	@Override
	public void begin() throws OutputException {
        if (!overwrite && outputFile.exists())
        	throw new OutputException("Output file already exists");
        try {
            outputFile.createNewFile();
        } catch (IOException e) {
        	throw new OutputException("Unable to create file for output", e);
        }
        if (!outputFile.canWrite())
        	throw new OutputException("Unable to write to output file");
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
    public void finalize() throws OutputException, OutputAlreadyWrittenException {
    	if (!outputWritten)
			end();
    }

    @Override
    public void end() throws OutputException, OutputAlreadyWrittenException {
    	if (outputWritten)
    		throw new OutputAlreadyWrittenException("The output has already been written to the file");
        List<MLArray> matlabData = MatlabDataExporter.buildMatlabSimulationList(plotDataListMap);
//        if (getRows() > 0)
//            matlabData.add(MatlabDataExporter.buildDouble("rows", getRows()));
//        if (getCols() > 0)
//            matlabData.add(MatlabDataExporter.buildDouble("cols", getCols()));
        try {
        	MatlabDataExporter.writeMatlabDataToFile(outputFile, matlabData);
            outputWritten = true;
        } catch (IOException e) {
        	throw new OutputException(e);
        }
    }

//    public int getRows() {
//        return rows;
//    }
//
//    public void setRows(int rows) {
//        this.rows = rows;
//    }
//
//    public int getCols() {
//        return cols;
//    }
//
//    public void setCols(int cols) {
//        this.cols = cols;
//    }

}
