package ch.ethz.khammash.ode;

public interface TimepointProvider {

	static class NoMoreTimepointsException extends RuntimeException {
		private static final long serialVersionUID = 4506876838123448120L;

		public NoMoreTimepointsException(String message) {
			super(message);
		}
	}

	double getInitialTimepoint();

	double getLastTimepoint();

	void reset();

	double getCurrentTimepoint();

	boolean hasNextTimepoint(double tCurrent);

	double getNextTimepoint(double tCurrent) throws NoMoreTimepointsException;

}
