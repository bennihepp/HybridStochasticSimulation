/**
 * @author      Benjamin Hepp <benjamin.hepp@bsse.ethz.ch>
 */
package ch.ethz.bhepp.ode;

/**
 * An Interface TimepointProvider to provide timesteps for an ODE {@link Solver}.
 */
public interface TimepointProvider {

	/**
	 * The Class NoMoreTimepointsException.
	 */
	static class NoMoreTimepointsException extends RuntimeException {
		
		/** The Constant serialVersionUID. */
		private static final long serialVersionUID = 4506876838123448120L;

		/**
		 * @param message the message
		 */
		public NoMoreTimepointsException(String message) {
			super(message);
		}
	}

	/**
	 * Gets the initial timepoint.
	 *
	 * @return the initial timepoint
	 */
	double getInitialTimepoint();

	/**
	 * Gets the last timepoint.
	 *
	 * @return the last timepoint
	 */
	double getLastTimepoint();

	/**
	 * Reset the provider.
	 */
	void reset();

	/**
	 * Gets the current timepoint.
	 *
	 * @return the current timepoint
	 */
	double getCurrentTimepoint();

	/**
	 * Checks if there is a timepoint > tCurrent.
	 *
	 * @param tCurrent the current time
	 * @return true, if there is a timepoint > tCurrent
	 */
	boolean hasNextTimepoint(double tCurrent);

	/**
	 * Gets the next timepoint > tCurrent.
	 *
	 * @param tCurrent the current time
	 * @return the next timepoint > tCurrent
	 * @throws NoMoreTimepointsException if there are no more timepoints > tCurrent
	 */
	double getNextTimepoint(double tCurrent) throws NoMoreTimepointsException;

}
