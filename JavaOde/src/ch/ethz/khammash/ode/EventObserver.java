/**
 * @author      Benjamin Hepp <benjamin.hepp@bsse.ethz.ch>
 */
package ch.ethz.khammash.ode;

/**
 * An Interface for observing the events/roots during the integration of an ODE with a {@link Solver}.
 */
public interface EventObserver {

	/**
	 * An Enum describing the Action to take when an event has been detected.
	 */
	enum EventAction {
		/** Stop the integration. */
		STOP, 
		/** Continue the integration. */
		CONTINUE,
	}

	/**
	 * Initialize the event observer (this is called by the ODE {@link Solver}, not the user).
	 *
	 * @param t0 the start time of integration
	 * @param x0 the initial state
	 * @param t1 the stop time of integration
	 */
	void initialize(double t0, double[] x0, double t1);

	/**
	 * Report an event that was encountered during ODE integration (this is called by the ODE {@link Solver}, not the user).
	 *
	 * @param eventIndex the index of the event
	 * @param t the time
	 * @param x the state
	 *
	 * @return the action to be performed
	 */
	EventAction report(int eventIndex, double t, double[] x);

}
