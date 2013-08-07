package ch.ethz.khammash.ode;

import java.nio.DoubleBuffer;


public class DirectBufferEventFunctionAdapter implements EventFunction {

	private EventFunction ef;
	private double[] x;
	private double[] values;
	private DoubleBuffer xBuffer;
	private DoubleBuffer valuesBuffer;

	public DirectBufferEventFunctionAdapter(EventFunction ef, int vectorFieldDimension,
			DoubleBuffer xBuffer, DoubleBuffer valuesBuffer) {
		this.ef = ef;
		this.x = new double[vectorFieldDimension];
		this.values = new double[ef.getNumberOfEventValues()];
		this.xBuffer = xBuffer;
		this.valuesBuffer = valuesBuffer;
	}

	@Override
	public int getNumberOfEventValues() {
    	return ef.getNumberOfEventValues();
    }

    public void computeEventValues(double t) {
		xBuffer.position(0);
		xBuffer.get(x);
		ef.computeEventValues(t, x, values);
		valuesBuffer.position(0);
		valuesBuffer.put(values);
    }

	@Override
    public void computeEventValues(double t, double[] x, double[] values) {
		ef.computeEventValues(t, x, values);
    }

}
