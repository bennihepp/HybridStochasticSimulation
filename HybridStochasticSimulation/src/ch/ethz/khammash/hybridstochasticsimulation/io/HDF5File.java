package ch.ethz.khammash.hybridstochasticsimulation.io;

import java.io.File;
import java.util.List;

import ncsa.hdf.object.Dataset;
import ncsa.hdf.object.Datatype;
import ncsa.hdf.object.FileFormat;
import ncsa.hdf.object.Group;
import ncsa.hdf.object.HObject;
import ncsa.hdf.object.h5.H5File;
import ncsa.hdf.object.h5.H5ScalarDS;

import com.google.common.base.Supplier;
import com.google.common.base.Suppliers;

public class HDF5File {

	public static HDF5File openForWriting(File file) throws HDF5Exception {
		return new HDF5File(openH5File(file));
	}

	public static HDF5File createForWriting(File file) throws HDF5Exception {
		return new HDF5File(createH5File(file));
	}

	public static class HDF5Exception extends Exception {

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

	public class HDF5Group {

		private final Group group;

		private HDF5Group(Group group) {
			this.group = group;
		}

		@Override
		public String toString() {
			return group.getFullName();
		}

		public HDF5Group createGroup(String name) throws HDF5Exception {
			try {
				return new HDF5Group(h5file.createGroup(name, this.group));
			} catch (Exception e) {
				throw new HDF5Exception(String.format("Unable to create child group for %s", this), e);
			}
		}

		public HDF5Group getGroup(String name) throws HDF5Exception {
			HObject member = getMember(name);
			if (member instanceof Group)
				return new HDF5Group((Group)member);
			throw new HDF5Exception(String.format("No such group: %s", name));
		}

		public HDF5Dataset getDataset(String name) throws HDF5Exception {
			HObject member = getMember(name);
			if (member instanceof Dataset)
				return new HDF5Dataset((Dataset)member);
			throw new HDF5Exception(String.format("No such dataset: %s", name));
		}

		public boolean hasGroup(String name) {
			HObject member = getMemberWithoutException(name);
			if (member instanceof Group)
				return true;
			return false;
		}

		public boolean hasDataset(String name) {
			HObject member = getMemberWithoutException(name);
			if (member instanceof Dataset)
				return true;
			return false;
		}

		private HObject getMember(String name) throws HDF5Exception {
			HObject member = getMemberWithoutException(name);
			if (member == null)
				throw new HDF5Exception(String.format("No such member: %s", name));
			return member;
		}

		private HObject getMemberWithoutException(String name) {
			List<HObject> members = group.getMemberList();
			for (HObject member : members) {
				if (member.getName().equals(name))
					return member;
			}
			return null;
		}

//		private boolean hasMember(String name) {
//			return getMemberWithoutException(name) != null;
//		}

		public Group getNativeGroup() {
			return group;
		}

		public HDF5Dataset createDataset(String name, HDF5Datatype dtype, HDF5Dimension dim) throws HDF5Exception {
			return createDataset(name, dtype, dim, null);
		}

		public HDF5Dataset createDataset(String name, HDF5Datatype dtype, HDF5Dimension dim, int gzip) throws HDF5Exception {
			return createDataset(name, dtype, dim, null, gzip);
		}

		public HDF5Dataset createDataset(String name, HDF5Datatype dtype, HDF5Dimension dim, double[] data) throws HDF5Exception {
			return createDataset(name, dtype, dim, data, 0);
		}

		public HDF5Dataset createDataset(String name, HDF5Datatype dtype, HDF5Dimension dim, double[] data, int gzip) throws HDF5Exception {
			try {
				Dataset dataset = h5file.createScalarDS(name, getNativeGroup(), dtype.getNativeDatatype(), dim.getDimensions(), dim.getMaxDimensions(), dim.getChunks(), gzip, data);
				return new HDF5Dataset(dataset);
			} catch (Exception e) {
				throw new HDF5Exception(String.format("Failed to create dataset %s", name), e);
			}
		}

//		public HDF5Dataset createGzipDataset(String name, HDF5Datatype dtype, HDF5Dimension dim, int gzipLevel) throws HDF5Exception {
//			checkArgument(0 <= gzipLevel && gzipLevel <= 9, "Expected 0 <= gzipLevel <= 9");
//			// Create the dataset creation property list, add the gzip compression filter
//	        try {
//	            int dcplId = H5.H5Pcreate(HDF5Constants.H5P_DATASET_CREATE);
//	            if (dcplId >= 0) {
//	                H5.H5Pset_deflate(dcplId, gzipLevel);
//	                // Set the chunk size
//	                H5.H5Pset_chunk(dcplId, dim.getNumberOfDimensions(), dim.getChunks());
//	            }
//
//	            // Create dataspace. Setting maximum size to NULL sets the maximum size to be the current size.
//                int filespaceId = H5.H5Screate_simple(dim.getNumberOfDimensions(), dim.getDimensions(), dim.getMaxDimensions());
//
//	            // Create the dataset
//	            int fileId = h5file.getFID();
//	            if ((fileId >= 0) && (filespaceId >= 0) && (dcplId >= 0)) {
//	                int datasetId = H5.H5Dcreate(
//	                		fileId, name, dtype.getNativeDatatype().toNative(), filespaceId,
//	                		HDF5Constants.H5P_DEFAULT, dcplId, HDF5Constants.H5P_DEFAULT);
//	            }
//	        catch (Exception e) {
//				throw new HDF5Exception(String.format("Failed to create gzip dataset %s", name));
//	        }
//		}

	}

	public class HDF5Datatype {

		private final Datatype dtype;

		private HDF5Datatype(Datatype dtype) {
			this.dtype = dtype;
		}

		@Override
		public String toString() {
			return dtype.getFullName();
		}

		public Datatype getNativeDatatype() {
			return dtype;
		}

	}

	public class HDF5Dataset {

		private final Dataset dataset;

		private HDF5Dataset(Dataset dataset) {
			this.dataset = dataset;
		}

		@Override
		public String toString() {
			return dataset.getFullName();
		}

		public void init() {
			dataset.init();
		}

		public Dataset getNativeDataset() {
			return dataset;
		}

		public long[] getDims() {
			return dataset.getDims();
		}

		public long[] getSelectedDims() {
			return dataset.getSelectedDims();
		}

		public long[] getStartDims() {
			return dataset.getStartDims();
		}

		public void extend(HDF5Dimension newDims) throws HDF5Exception {
			try {
				((H5ScalarDS)dataset).extend(newDims.getDimensions());
			} catch (ncsa.hdf.hdf5lib.exceptions.HDF5Exception e) {
				throw new HDF5Exception("Unable to extend dataset", e);
			}
		}

		public void write(Object data) throws HDF5Exception {
			try {
				dataset.write(data);
			} catch (Exception e) {
				throw new HDF5Exception("Unable to write to dataset", e);
			}
		}

	}

	private static class HDF5FileFormat {

		private final FileFormat ff;

		public static HDF5FileFormat createHDF5FileFormat() throws HDF5Exception {
	    	FileFormat ff = FileFormat.getFileFormat(FileFormat.FILE_TYPE_HDF5);
	    	if (ff == null)
	    		throw new HDF5Exception("Cannot find HDF5 FileFormat");

	    	return new HDF5FileFormat(ff);
		}

		public HDF5FileFormat(FileFormat ff) {
			this.ff = ff;
		}

		public H5File openNativeFile(File file, boolean truncate) throws HDF5Exception {
	    	try {
	    		H5File h5File;
	    		// The default HDF5 file creation and access properties are used.
	    		if (truncate)
		    		// If the specified file already exists, it is truncated.
	    			h5File = (H5File)ff.createFile(file.getAbsolutePath(), FileFormat.FILE_CREATE_DELETE);
	    		else
		    		h5File = (H5File)ff.createFile(file.getAbsolutePath(), FileFormat.FILE_CREATE_OPEN);

	    		// Check for error condition and report.
	    		if (h5File == null)
	    			throw new HDF5Exception(String.format("Failed to create file: {}", file.getAbsolutePath()));

	    		return h5File;

	    	} catch (Exception e) {
	    		throw new HDF5Exception("Unable to create HDF5 file", e);
	    	}
		}

	}

	private class RootGroupSupplier implements Supplier<HDF5Group> {

		@Override
		public HDF5Group get() {
			Group rootGroup = (Group)((javax.swing.tree.DefaultMutableTreeNode)h5file.getRootNode()).getUserObject();
			return new HDF5Group(rootGroup);
		}

	}

	private final H5File h5file;
	private final Supplier<HDF5Group> memoizedRootGroup;

	private HDF5File(H5File h5file) {
		this.h5file = h5file;
		this.memoizedRootGroup = Suppliers.memoize(new RootGroupSupplier());
	}

	public void open() throws HDF5Exception {
		try {
			h5file.open();
		} catch (Exception e) {
			throw new HDF5Exception("Unable to open HDF5 file", e);
		}
	}

	public HDF5Group getRootGroup() {
		return memoizedRootGroup.get();
	}

    private static H5File openH5File(File file) throws HDF5Exception {
    	HDF5FileFormat fileFormat = HDF5FileFormat.createHDF5FileFormat();
    	return fileFormat.openNativeFile(file, false);
	}

    private static H5File createH5File(File file) throws HDF5Exception {
    	HDF5FileFormat fileFormat = HDF5FileFormat.createHDF5FileFormat();
    	return fileFormat.openNativeFile(file, true);
	}

//    private static H5File createGzipH5File(File file) throws HDF5Exception {
//    	if (!checkGzipFilterAvailable())
//    		throw new HDF5Exception("Gzip filter is not available");
//    	HDF5FileFormat fileFormat = HDF5FileFormat.createHDF5FileFormat();
//    	return fileFormat.createNativeFile(file);
//	}

//	private static boolean checkGzipFilterAvailable() throws HDF5Exception {
//        try {
//            int available = H5.H5Zfilter_avail(HDF5Constants.H5Z_FILTER_DEFLATE);
//            if (available == 0)
//                return false;
//
//            int filter_info = H5.H5Zget_filter_info(HDF5Constants.H5Z_FILTER_DEFLATE);
//            if (((filter_info & HDF5Constants.H5Z_FILTER_CONFIG_ENCODE_ENABLED) == 0)
//                    || ((filter_info & HDF5Constants.H5Z_FILTER_CONFIG_DECODE_ENABLED) == 0)) {
//            	return false;
//            }
//        }
//        catch (Exception e) {
//        	throw new HDF5Exception("Error while trying to check for gzip filter");
//        }
//        return true;
//	}

	public HDF5Datatype getDoubleDatatype() throws HDF5Exception {
		return getDatatype(Datatype.CLASS_FLOAT, 8, Datatype.NATIVE, -1);
	}

	public HDF5Datatype getDatatype(int tclass, int tsize, int tsign, int torder) throws HDF5Exception {
		Datatype dtype;
		try {
			dtype = h5file.createDatatype(tclass, tsize, torder, tsign);
		} catch (Exception e) {
			String msg = String.format(
					"Unable to create Datatype (tclass=%d, tsize=%d, torder=%d, tsign=%d",
					tclass, tsize, torder, tsign);
			throw new HDF5Exception(msg, e);
		}
		return new HDF5Datatype(dtype);
	}

	public void close() throws HDF5Exception {
		try {
			h5file.close();
		} catch (ncsa.hdf.hdf5lib.exceptions.HDF5Exception e) {
			throw new HDF5Exception("Unable to close HDF5 file", e);
		}
	}

}
