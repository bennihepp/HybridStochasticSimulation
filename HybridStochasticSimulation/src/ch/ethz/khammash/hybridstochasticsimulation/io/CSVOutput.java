package ch.ethz.khammash.hybridstochasticsimulation.batch;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import au.com.bytecode.opencsv.CSVWriter;
import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob.OutputAlreadyWrittenException;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FinitePlotData;


public class CSVOutput implements SimulationOutput {

    private File outputFile;
	private boolean overwrite;
	private CSVWriter writer;
    private boolean outputWritten = false;

    public CSVOutput(String outputFilename, boolean overwrite) throws IOException {
        this(new File(outputFilename), overwrite);
    }

    public CSVOutput(File outputFile, boolean overwrite) throws IOException {
        this.outputFile = outputFile;
        this.overwrite = overwrite;
    }

    @Override
    public String toString() {
        return outputFile.getAbsolutePath();
    }

	@Override
	public void init() throws OutputException {
        if (!overwrite && outputFile.exists())
        	throw new OutputException("Output file already exists");
        try {
            outputFile.createNewFile();
        } catch (IOException e) {
        	throw new OutputException("Unable to create file for output", e);
        }
        if (!outputFile.canWrite())
        	throw new OutputException("Unable to write to output file");
        try {
        	writer = new CSVWriter(new FileWriter(outputFile));
        } catch (IOException e) {
        	throw new OutputException(e);
        }
	}

    @Override
    public void add(String simulationName, FinitePlotData plotData) {
    	writePlotData(simulationName, plotData);
    }

    @Override
    public void addAll(String simulationName, List<FinitePlotData> plotDataList) {
    	for (FinitePlotData plotData : plotDataList)
    		add(simulationName, plotData);
    }

    private void writePlotData(String simulationName, FinitePlotData plotData) {
    	String[] simulationColumns = {"Simulation", simulationName};
    	writer.writeNext(simulationColumns);
    	String[] sizeColumns = {"NumberOfStates", Integer.toString(plotData.getNumberOfStates()),
    			"NumberOfTimepoitns", Integer.toString(plotData.getNumberOfTimePoints())};
    	writer.writeNext(sizeColumns);
		String[] columns = new String[plotData.getNumberOfStates()];
    	for (int state=0; state < plotData.getNumberOfStates(); state++) {
    		for (int i=0; i < plotData.getNumberOfStates(); i++)
    			columns[i] = Double.toString(plotData.getxState(state, i));
    		writer.writeNext(columns);
    	}
    }

    @Override
    public void finalize() throws OutputException, OutputAlreadyWrittenException {
    	if (!outputWritten)
			write();
    }

    @Override
    public void write() throws OutputException, OutputAlreadyWrittenException {
    	if (outputWritten)
    		throw new OutputAlreadyWrittenException("The output has already been written to the file");
    	try {
			writer.close();
		} catch (IOException e) {
			throw new OutputException(e);
		}
    }

}
