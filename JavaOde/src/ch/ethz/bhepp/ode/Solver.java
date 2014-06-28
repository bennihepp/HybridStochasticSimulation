/**
 * @author      Benjamin Hepp <benjamin.hepp@bsse.ethz.ch>
 */
package ch.ethz.bhepp.ode;


/**
 * An Interface for solving an {@link Ode}.
 */
public interface Solver {

	/**
	 * The Class InitializationException.
	 */
	static class InitializationException extends RuntimeException {
		
		/** The Constant serialVersionUID. */
		private static final long serialVersionUID = 4775536474978013594L;

		/**
		 * @param message the message
		 */
		public InitializationException(String message) {
			super(message);
		}

		/**
		 * @param message the message
		 * @param cause the underlying cause of the exception
		 */
		public InitializationException(String message, Throwable cause) {
			super(message, cause);
		}
	}

	/**
	 * The Class IntegrationException.
	 */
	static class IntegrationException extends RuntimeException {
		
		/** The Constant serialVersionUID. */
		private static final long serialVersionUID = 8509338706972523801L;

		/**
		 * @param message the message
		 */
		public IntegrationException(String message) {
			super(message);
		}

		/**
		 * @param message the message
		 * @param cause the underlying cause of the exception
		 */
		public IntegrationException(String message, Throwable cause) {
			super(message, cause);
		}
	}

    /**
     * Initialize the solver.
     *
     * @param ode the ODE to solve
     * @param ef the function defining the roots to be found
     * @param stateObserver the observer of the states at the integration time-steps
     * @param eventObserver the observer of the root functions
     * @throws InitializationException if there was a problem during initialization
     */
    void initialize(Ode ode, EventFunction ef, StateObserver stateObserver, EventObserver eventObserver) throws InitializationException;

    /**
     * Dispose the solver. This is important for solvers that make use of JNI.
     */
    void dispose();

    /**
     * Prepare a call to {@link #integrate()}.
     *
     * @param timepointProvider defines the time-steps for the integration
     * @param x0 the initial state
     */
    public void prepare(TimepointProvider timepointProvider, double[] x0);

    /**
     * Integrate the ODE. The user is responsible to call {@link #prepare} before calling this function.
     *
     * @return the time where the integration has stopped
     * @throws IntegrationException if there was a problem during the integration
     * @throws IllegalStateException if {@link #initialize} has not been called before
     */
    public double integrate() throws IntegrationException, IllegalStateException;

    /**
     * Integrate the ODE.
     *
     * @param t0 the start time of the integration
     * @param x0 the initial state
     * @param t1 the stop time of the integration
     * @return the time where the integration has stopped
     * @throws IntegrationException if there was a problem during the integration
     * @throws UnsupportedOperationException if the solver only implements integration with {@link TimepointProvider}s
     * @throws IllegalStateException if {@link #initialize} has not been called before
     */
    public double integrate(double t0, double[] x0, double t1) throws IntegrationException, UnsupportedOperationException, IllegalStateException;

    /**
     * Integrate the ODE.
     *
     * @param timepointProvider defines the time-steps for the integration
     * @param x0 the initial state
     * @return the time where the integration has stopped
     * @throws IntegrationException if there was a problem during the integration
     * @throws IllegalStateException if {@link #initialize} has not been called before
     */
    public double integrate(TimepointProvider timepointProvider, double[] x0) throws IntegrationException, IllegalStateException;

}
