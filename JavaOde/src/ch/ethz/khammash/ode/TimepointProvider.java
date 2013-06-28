package ch.ethz.khammash.ode;

public interface TimepointProvider {

	public double getInitialTimepoint();

	public double getLastTimepoint();

	public void reset();

	public double getCurrentTimepoint();

	public boolean hasNextTimepoint();

	public double getNextTimepoint();

}
