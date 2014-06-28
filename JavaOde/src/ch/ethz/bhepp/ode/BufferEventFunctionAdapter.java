/**
 * @author      Benjamin Hepp <benjamin.hepp@bsse.ethz.ch>
 */
package ch.ethz.bhepp.ode;

import java.nio.DoubleBuffer;


/**
 * This class takes an {@link EventFunction} and provides an equivalent {@link EventFunction} implementation that is backed by
 * buffers suitable for solvers that make use of JNI (by using direct buffers).
 */
public class BufferEventFunctionAdapter implements EventFunction {

	/** The {@link EventFunction}. */
	private EventFunction ef;

	/** An array for the state. */
	private double[] x;

	/** An array for the event/root values. */
	private double[] values;

	/** A buffer for the state. */
	private DoubleBuffer xBuffer;

	/** A buffer for the event/root values. */
	private DoubleBuffer valuesBuffer;

	/**
	 * Instantiates a new direct buffer adapter for an {@link EventFunction}.
	 *
	 * @param ef the {@link EventFunction}
	 * @param vectorFieldDimension the dimension of the {@link Ode} vector field
	 * @param xBuffer the buffer for the state
	 * @param valuesBuffer the buffer for the event/root values (will be modified by later calls)
	 */
	public BufferEventFunctionAdapter(EventFunction ef, int vectorFieldDimension,
			DoubleBuffer xBuffer, DoubleBuffer valuesBuffer) {
		this.ef = ef;
		this.x = new double[vectorFieldDimension];
		this.values = new double[ef.getNumberOfEventValues()];
		this.xBuffer = xBuffer;
		this.valuesBuffer = valuesBuffer;
	}

	/* (non-Javadoc)
	 * @see ch.ethz.bhepp.ode.EventFunction#getNumberOfEventValues()
	 */
	@Override
	public int getNumberOfEventValues() {
    	return ef.getNumberOfEventValues();
    }

    /**
	 * Compute the event/root values based on the time t and the state in the state buffer
	 * and store it in the event/root values buffer
	 *
	 * @param t the time
     */
    public void computeEventValues(double t) {
		xBuffer.position(0);
		xBuffer.get(x);
		ef.computeEventValues(t, x, values);
		valuesBuffer.position(0);
		valuesBuffer.put(values);
    }

	/* (non-Javadoc)
	 * @see ch.ethz.bhepp.ode.EventFunction#computeEventValues(double, double[], double[])
	 */
	@Override
    public void computeEventValues(double t, double[] x, double[] values) {
		ef.computeEventValues(t, x, values);
    }

}
