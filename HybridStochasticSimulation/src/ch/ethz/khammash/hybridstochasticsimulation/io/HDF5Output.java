package ch.ethz.khammash.hybridstochasticsimulation.io;

import static com.google.common.base.Preconditions.checkArgument;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob.OutputAlreadyWrittenException;
import ch.ethz.khammash.hybridstochasticsimulation.io.HDF5File.HDF5Datatype;
import ch.ethz.khammash.hybridstochasticsimulation.io.HDF5File.HDF5Group;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FinitePlotData;


// TODO: Add metadata to HDF5 file (number of simulations, species labels, etc.)
public class HDF5Output implements SimulationOutput {

	private static final Logger logger = LoggerFactory.getLogger(HDF5Output.class);

    private File outputFile;
	private boolean overwrite;
    private boolean outputWritten = false;
	private int chunkSize = 64;
	private int gzipLevel = 0;
	private HDF5File h5file;
//	private H5File h5file;
	private HDF5Group simulationsGroup;
	private Map<String, HDF5File.HDF5Dataset> datasetMap;

    public HDF5Output(String outputFilename, boolean overwrite) throws IOException {
        this(new File(outputFilename), overwrite);
    }

    public HDF5Output(File outputFile, boolean overwrite) throws IOException {
        this.outputFile = outputFile;
        this.overwrite = overwrite;
        datasetMap = new HashMap<>();
    }

    @Override
    public String toString() {
        return outputFile.getAbsolutePath();
    }

    public void setChunkSize(int chunkSize) {
    	checkArgument(1 <= chunkSize, "chunk size has to be greater equal than 1");
    	this.chunkSize = chunkSize;
    }

    public void setGzipLevel(int gzipLevel) {
    	checkArgument(0 <= gzipLevel && gzipLevel <= 9,
    			"gzip level has to be between 0 and 9");
    	this.gzipLevel = gzipLevel;
    }

	@Override
	public void begin() throws OutputException {
    	if (outputWritten)
    		throw new OutputAlreadyWrittenException("The output has already been written to the file");
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

    		h5file = HDF5File.openForWriting(outputFile);
    		h5file.open();

    		 HDF5Group rootGroup = h5file.getRootGroup();
    		 simulationsGroup = rootGroup.createGroup("simulations");

		} catch (Exception e) {
			throw new OutputException(e);
		}
	}

//    private H5File createHDF5File(File outputFile) throws Exception {
//    	// Retrieve an instance of the implementing class for the HDF5 format
//    	FileFormat fileFormat = FileFormat.getFileFormat(FileFormat.FILE_TYPE_HDF5);
//    	if (fileFormat == null)
//    		throw new HDF5Exception("Cannot find HDF5 FileFormat");
//
//		// If the specified file already exists, it is truncated.
//		// The default HDF5 file creation and access properties are used.
//		H5File file = (H5File)fileFormat.createFile(outputFile.getAbsolutePath(), FileFormat.FILE_CREATE_DELETE);
//		
//		// Check for error condition and report.
//		if (file == null)
//			throw new HDF5Exception(String.format("Failed to create file: {}", outputFile.getAbsolutePath()));
//		
//		return file;
//	}

	@Override
    public void add(String simulationName, FinitePlotData plotData) throws OutputException {
    	writePlotData(simulationName, plotData);
    }

    @Override
    public void addAll(String simulationName, List<FinitePlotData> plotDataList) throws OutputException {
    	for (FinitePlotData plotData : plotDataList)
    		add(simulationName, plotData);
    }

    private void writePlotData(String simulationName, FinitePlotData plotData) throws OutputException {
    	try {
    		HDF5File.HDF5Dataset xSeriesDataset;
    		if (datasetMap.containsKey(simulationName)) {
    			xSeriesDataset = datasetMap.get(simulationName);
    		} else {
    			HDF5Group simulationGroup = simulationsGroup.createGroup(simulationName);
    			HDF5Datatype dtype = h5file.getDoubleDatatype();
				double[] tSeries = plotData.gettSeries();
				long[] tSeriesDims1D = { tSeries.length };
				HDF5Dimension tSeriesDim = new HDF5Dimension(tSeriesDims1D);
				simulationGroup.createDataset("tSeries", dtype, tSeriesDim, tSeries);
				long[] xSeriesDims3D = { 0, plotData.getNumberOfStates(), plotData.getNumberOfTimePoints() };
				long[] xSeriesMaxDims3D = { Integer.MAX_VALUE, plotData.getNumberOfStates(), plotData.getNumberOfTimePoints() };
				long[] xSeriesChunks = xSeriesDims3D.clone();
				xSeriesChunks[0] = chunkSize;
				HDF5Dimension xSeriesDim = new HDF5Dimension(xSeriesDims3D, xSeriesMaxDims3D, xSeriesChunks);
    	    	if (logger.isDebugEnabled()) {
    	    		logger.debug("Creating dataset with dimension={}, gzipLevel={}", xSeriesDim, gzipLevel);
    	    	}
    	    	xSeriesDataset = simulationGroup.createDataset("xSeries", dtype, xSeriesDim, gzipLevel);
//				xSeriesDataset = (H5ScalarDS)h5file.createScalarDS("xSeries", simulationGroup, dtype, xSeriesDims3D, xSeriesMaxDims3D, xSeriesChunks, 0, null);
				xSeriesDataset.init();
				long[] selected = xSeriesDataset.getSelectedDims();
				selected[0] = 1;
				selected[1] = xSeriesDims3D[1];
				selected[2] = xSeriesDims3D[2];
				long[] start = xSeriesDataset.getStartDims();
				start[0]--;
				datasetMap.put(simulationName, xSeriesDataset);
    		}
			// extend dataset
			long[] xSeriesDims3D = xSeriesDataset.getDims();
			long[] newDims3D = xSeriesDims3D.clone();
			newDims3D[0]++;
	    	if (logger.isDebugEnabled()) {
	    		long x = newDims3D[0];
	    		long y = newDims3D[1];
	    		long z = newDims3D[2];
	    		logger.debug("Extending dataset to size=({},{},{})", x, y, z);
	    	}
	    	HDF5Dimension newDims = new HDF5Dimension(newDims3D);
			xSeriesDataset.extend(newDims);
			long[] start = xSeriesDataset.getStartDims();
			start[0]++;
			// TOOD: probably data can be used directly here
	        double[][][] data = new double[1][][];
	        data[0] = plotData.getxSeries();
	    	if (logger.isDebugEnabled())
	    		logger.debug("Writing trajectory to file");
	        xSeriesDataset.write(data);
		} catch (Exception e) {
			throw new OutputException(e);
		}
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
    	try {
			h5file.close();
		} catch (Exception e) {
			throw new OutputException("Error while closing file", e);
		}
    }

}
