package ch.ethz.khammash.hybridstochasticsimulation.io;

import static com.google.common.base.Preconditions.checkArgument;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ncsa.hdf.object.Datatype;
import ncsa.hdf.object.FileFormat;
import ncsa.hdf.object.Group;
import ncsa.hdf.object.h5.H5File;
import ncsa.hdf.object.h5.H5ScalarDS;
import ch.ethz.khammash.hybridstochasticsimulation.batch.SimulationJob.OutputAlreadyWrittenException;
import ch.ethz.khammash.hybridstochasticsimulation.trajectories.FinitePlotData;


public class HDF5Output implements SimulationOutput {

	public static class HDF5Exception extends OutputException {

		private static final long serialVersionUID = 1L;

		public HDF5Exception(String msg) {
			super(msg);
		}

		public HDF5Exception(String msg, Throwable cause) {
			super(msg, cause);
		}

		public HDF5Exception(Throwable cause) {
			super(cause);
		}

	}

    private File outputFile;
	private boolean overwrite;
    private boolean outputWritten = false;
	private int gzipLevel = 0;
	private H5File h5file;
	private Group simulationsGroup;
	private Map<String, H5ScalarDS> datasetMap;

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

    		h5file = createHDF5File(outputFile);
			h5file.open();

			Group rootGroup = (Group)((javax.swing.tree.DefaultMutableTreeNode)h5file.getRootNode()).getUserObject();
			simulationsGroup = h5file.createGroup("simulations", rootGroup);

		} catch (Exception e) {
			throw new HDF5Exception(e);
		}
	}

    private H5File createHDF5File(File outputFile) throws Exception {
    	// Retrieve an instance of the implementing class for the HDF5 format
    	FileFormat fileFormat = FileFormat.getFileFormat(FileFormat.FILE_TYPE_HDF5);
    	if (fileFormat == null)
    		throw new HDF5Exception("Cannot find HDF5 FileFormat");

     // If the specified file already exists, it is truncated.
     // The default HDF5 file creation and access properties are used.
     H5File file = (H5File)fileFormat.createFile(outputFile.getAbsolutePath(), FileFormat.FILE_CREATE_DELETE);

     // Check for error condition and report.
     if (file == null)
    	 throw new HDF5Exception(String.format("Failed to create file: {}", outputFile.getAbsolutePath()));

     return file;
	}

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
    		H5ScalarDS xSeriesDataset;
    		if (datasetMap.containsKey(simulationName)) {
    			xSeriesDataset = datasetMap.get(simulationName);
    		} else {
    			Group simulationGroup = h5file.createGroup(simulationName, simulationsGroup);
				Datatype dtype = h5file.createDatatype(Datatype.CLASS_FLOAT, 8, Datatype.NATIVE, -1);
				double[] tSeries = plotData.gettSeries();
				long[] tSeriesDims1D = { tSeries.length };
				h5file.createScalarDS("tSeries", simulationGroup, dtype, tSeriesDims1D, null, null, gzipLevel, tSeries);
				long[] xSeriesDims3D = { 0, plotData.getNumberOfStates(), plotData.getNumberOfTimePoints() };
				long[] xSeriesMaxDims3D = { Integer.MAX_VALUE, plotData.getNumberOfStates(), plotData.getNumberOfTimePoints() };
				long[] xSeriesChunks = xSeriesDims3D.clone();
				xSeriesChunks[0]++;
				xSeriesDataset = (H5ScalarDS)h5file.createScalarDS("xSeries", simulationGroup, dtype, xSeriesDims3D, xSeriesMaxDims3D, xSeriesChunks, 0, null);
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
			xSeriesDataset.extend(newDims3D);
			long[] start = xSeriesDataset.getStartDims();
			start[0]++;
			// TOOD: probably data can be used directly here
	        double[][][] data = new double[1][][];
	        data[0] = plotData.getxSeries();
	        xSeriesDataset.write(data);
		} catch (Exception e) {
			throw new HDF5Exception(e);
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
			throw new HDF5Exception("Error while closing file", e);
		}
    }

}
