/**
 * @author      Benjamin Hepp <benjamin.hepp@bsse.ethz.ch>
 */
package ch.ethz.bhepp.ode;

import ch.ethz.bhepp.ode.Solver.InitializationException;


/**
 * An Interface for solving an {@link Ode} with adaptive step size.
 */
public interface AdaptiveStepSolver {

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
     * @param xOut array to store the new state after the integration step
     * 
     * @return the new timepoint after the integration step
     *
     */
	double integrateStep(double[] xOut);

    /**
     * Compute the interpolated solution at a timepoint within the last integration step. {@link #integrateStep} has to be called before calling this method.
     * 
     * @param t the timepoint at which the interpolated solution should be computed
     * @param xOut array to store the new state after the integration step
     *
     * @throws IllegalStateException if {@link #integrateStep} has not been called before, the last integration step failed or the solution could not be interpolated.
     */
	void computeInterpolatedSolution(double t, double[] xOut) throws IllegalStateException;

    /**
     * Get the current step size of the solver (used for the next step). {@link #integrateStep} has to be called before calling this method.
     * 
     * @throws IllegalStateException if {@link #initialize} has not been called before
     * @throws UnsupportedOperationException if the solver doesn't support retrieving the current step size
     */
	double getCurrentStepSize() throws IllegalStateException, UnsupportedOperationException;

    /**
     * Get the last step size of the solver (used for the last step). {@link #integrateStep} has to be called before calling this method.
     * 
     * @throws IllegalStateException if {@link #initialize} has not been called before
     * @throws UnsupportedOperationException if the solver doesn't support retrieving the last step size
     */
	double getLastStepSize() throws IllegalStateException, UnsupportedOperationException;

    /**
     * Set the current step size for the next call to {@link #integrateStep}.
     * 
     * @throws IllegalStateException if {@link #initialize} has not been called before
     * @throws UnsupportedOperationException if the solver doesn't support modifying the step size
     */
	void setCurrentStepSize(double stepSize) throws IllegalStateException, UnsupportedOperationException;

    /**
     * Set the minimum step size for the next call to {@link #integrateStep}.
     * 
     * @throws IllegalStateException if {@link #initialize} has not been called before
     * @throws UnsupportedOperationException if the solver doesn't support modifying the minimum step size
     */
	void setMinimumStepSize(double minStepSize) throws IllegalStateException, UnsupportedOperationException;

    /**
     * Set the maximum step size for the next call to {@link #integrateStep}.
     * 
     * @throws IllegalStateException if {@link #initialize} has not been called before
     * @throws UnsupportedOperationException if the solver doesn't support modifying the maximum step size
     */
	void setMaximumStepSize(double maxStepSize) throws IllegalStateException, UnsupportedOperationException;

}
