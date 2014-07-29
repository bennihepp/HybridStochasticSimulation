/**
 * @author      Benjamin Hepp <benjamin.hepp@bsse.ethz.ch>
 */
package ch.ethz.bhepp.ode;

import ch.ethz.bhepp.ode.Solver.InitializationException;


/**
 * An Interface for solving an {@link Ode} with fixed step size.
 */
public interface FixedStepSolver {

    /**
     * Initialize the solver.
     *
     * @param ode the ODE to solve
     * @throws InitializationException if there was a problem during initialization
     */
    void initialize(Ode ode) throws InitializationException;

    /**
     * Prepare a call to {@link #integrateStep()}.
     *
     * @param t0 defines the current timepoint for the integration
     * @param x0 defines the current state for the integration
     * @param t1 defines the end timepoint for the integration
     * 
     * @throws IllegalStateException if {@link #initialize} has not been called before
     */
	void prepareStep(double t0, double[] x0, double t1) throws IllegalStateException;

    /**
     * Integrate the ODE with a single step. {@link #prepareStep} has to be called before calling this method.
     * 
     * @return the new timepoint after the integration step
     *
     */
	double integrateStep();

    /**
     * Integrate the ODE with a single step. {@link #prepareStep} has to be called before calling this method.
     *
     * @param step the size of the integration step to take
     * 
     * @return the new timepoint after the integration step
     */
	double integrateStep(double step);

    /**
     * Integrate the ODE with a single step.
     *
     * @param t current timepoint for the integration
     * @param x array holding the current state for the integration
     * @param t1 defines the end timepoint for the integration
     * @param xOut array to store the new state after the integration step
     * 
     * @return the new timepoint after the integration step
     */
	double integrateStep(double t, double[] x, double t1, double[] xOut);

    /**
     * Integrate the ODE with a single step.
     *
     * @param t current timepoint for the integration
     * @param x array holding the current state for the integration
     * @param xOut array to store the new state after the integration step
     * @param step the size of the integration step to take
     * 
     * @return the new timepoint after the integration step
     */
	double integrateStep(double t, double[] x, double[] xOut, double step);

	/**
	 * Set step size of integration.
	 * 
	 * @param stepSize new step size for integration
	 * @throws UnsupportedOperationException if changing step size is not possible
	 */
	void setStepSize(double stepSize) throws UnsupportedOperationException;
}
