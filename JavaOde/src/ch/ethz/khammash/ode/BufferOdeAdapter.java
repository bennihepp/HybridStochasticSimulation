/**
 * @author      Benjamin Hepp <benjamin.hepp@bsse.ethz.ch>
 */
package ch.ethz.khammash.ode;

import java.nio.DoubleBuffer;


/**
 * This class takes an {@link Ode} and provides an equivalent {@link Ode} implementation that is backed by
 * buffers suitable for solvers that make use of JNI (by using direct buffers).
 */
public class BufferOdeAdapter implements Ode {

	/** The {@link Ode}. */
	private Ode ode;
	
	/** An array for the state. */
	private double[] x;
	
	/** An array for the vector field. */
	private double[] xDot;
	
	/** A buffer for the state. */
	private DoubleBuffer xBuffer;
	
	/** A buffer for the vector field. */
	private DoubleBuffer xDotBuffer;

	/**
	 * Instantiates a new direct buffer adapter for an {@link Ode}.
	 *
	 * @param ode the {@link Ode}
	 * @param xBuffer the buffer for the state
	 * @param xDotBuffer the buffer for the vector field (will be modified by later calls)
	 */
	public BufferOdeAdapter(Ode ode, DoubleBuffer xBuffer, DoubleBuffer xDotBuffer) {
		this.ode = ode;
		this.x = new double[ode.getDimensionOfVectorField()];
		this.xDot = new double[ode.getDimensionOfVectorField()];
		this.xBuffer = xBuffer;
		this.xDotBuffer = xDotBuffer;
	}

	/* (non-Javadoc)
	 * @see ch.ethz.khammash.ode.Ode#getDimensionOfVectorField()
	 */
	@Override
	public int getDimensionOfVectorField() {
		return ode.getDimensionOfVectorField();
	}

	/**
	 * Compute the vector field based on the time t and the state in the state buffer
	 * and store it in the vector field buffer.
	 *
	 * @param t the time
	 */
	public void computeVectorField(double t) {
		xBuffer.position(0);
		xBuffer.get(x);
		ode.computeVectorField(t, x, xDot);
		xDotBuffer.position(0);
		xDotBuffer.put(xDot);
	}

	/* (non-Javadoc)
	 * @see ch.ethz.khammash.ode.Ode#computeVectorField(double, double[], double[])
	 */
	@Override
	public void computeVectorField(double t, double[] x, double[] xDot) {
		ode.computeVectorField(t, x, xDot);
	}

}
