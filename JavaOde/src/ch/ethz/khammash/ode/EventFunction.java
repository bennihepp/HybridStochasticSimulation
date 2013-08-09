/**
 * @author      Benjamin Hepp <benjamin.hepp@bsse.ethz.ch>
 */
package ch.ethz.khammash.ode;

/**
 * An Interface defining the roots to be found during integration.
 */
public interface EventFunction {

    /**
     * Gets the number of event/root values.
     *
     * @return the number of event values
     */
    int getNumberOfEventValues();

    /**
     * Compute the event/root values.
     *
     * @param t the time
     * @param x the state (won't be modified!)
     * @param values the array for storing the event/root values (will be modified!)
     */
    void computeEventValues(double t, double[] x, double[] values);

}
