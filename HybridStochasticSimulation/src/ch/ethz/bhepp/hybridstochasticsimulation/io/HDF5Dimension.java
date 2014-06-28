package ch.ethz.bhepp.hybridstochasticsimulation.io;

import java.util.Arrays;

public class HDF5Dimension {

	private final long[] dim;
	private final long[] maxDim;
	private final long[] chunks;

	public HDF5Dimension(long[] dim) {
		this(dim, null, null);
	}

	public HDF5Dimension(long[] dim, long[] maxDim) {
		this(dim, maxDim, null);
	}

	public HDF5Dimension(long[] dim, long[] maxDim, long[] chunk) {
		this.dim = dim;
		this.maxDim = maxDim;
		this.chunks = chunk;
	}

	@Override
	public String toString() {
		return String.format("{%s, %s, %s}",
				Arrays.toString(dim),
				Arrays.toString(maxDim),
				Arrays.toString(chunks));
	}

	public long[] getDimensions() {
		return cloneIfNotNull(dim);
	}

	public long[] getMaxDimensions() {
		return cloneIfNotNull(maxDim);
	}

	public long[] getChunks() {
		return cloneIfNotNull(chunks);
	}

	public int getNumberOfDimensions() {
		return dim == null ? -1 : dim.length;
	}

	private long[] cloneIfNotNull(long[] arr) {
		if (arr == null)
			return null;
		return arr.clone();
	}

}
