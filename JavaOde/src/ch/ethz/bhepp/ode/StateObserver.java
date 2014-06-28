/**
 * @author      Benjamin Hepp <benjamin.hepp@bsse.ethz.ch>
 */
package ch.ethz.bhepp.ode;

/**
 * An Interface for observing the state during the integration of an ODE with a {@link Solver}.
 */
public interface StateObserver {

	/**
	 * Initialize the state observer (this is called by the ODE {@link Solver}, not the user).
	 *
	 * @param t0 the start time of integration
	 * @param x0 the initial state
	 * @param t1 the stop time of integration
	 */
	void initialize(double t0, double[] x0, double t1);

	/**
	 * Report a state that was encountered during ODE integration (this is called by the ODE {@link Solver}, not the user).
	 *
	 * @param t the time
	 * @param x the state
	 */
	void report(double t, double[] x);

}
