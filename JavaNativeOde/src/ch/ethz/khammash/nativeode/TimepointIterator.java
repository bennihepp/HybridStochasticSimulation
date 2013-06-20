package ch.ethz.khammash.nativeode;

public interface TimepointIterator {

	public double getInitialTimepoint();

	public double getLastTimepoint();

	public void reset();

	public double getCurrentTimepoint();

	public boolean hasNextTimepoint();

	public double getNextTimepoint();

}
