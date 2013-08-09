/**
 * @author      Benjamin Hepp <benjamin.hepp@bsse.ethz.ch>
 */
package ch.ethz.khammash.ode;

/**
 * An Interface describing a system of Ordinary Differential Equations.
 */
public interface Ode {

    /**
     * Gets the dimension of the vector field.
     *
     * @return the dimension of the vector field
     */
    int getDimensionOfVectorField();
    
    /**
     * Compute the vector field of the ODE (the time derivative).
     *
     * @param t the time
     * @param x the state (won't be modified!)
     * @param xDot the array for storing the vector field (will be modified!)
     */
    void computeVectorField(double t, double[] x, double[] xDot);

};
