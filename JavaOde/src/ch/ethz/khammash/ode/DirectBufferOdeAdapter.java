package ch.ethz.khammash.ode;

import java.nio.DoubleBuffer;


public class DirectBufferOdeAdapter implements Ode {

	private Ode ode;
	private double[] x;
	private double[] xDot;
	private DoubleBuffer xBuffer;
	private DoubleBuffer xDotBuffer;

	public DirectBufferOdeAdapter(Ode ode, DoubleBuffer xBuffer, DoubleBuffer xDotBuffer) {
		this.ode = ode;
		this.x = new double[ode.getDimensionOfVectorField()];
		this.xDot = new double[ode.getDimensionOfVectorField()];
		this.xBuffer = xBuffer;
		this.xDotBuffer = xDotBuffer;
	}

	@Override
	public int getDimensionOfVectorField() {
		return ode.getDimensionOfVectorField();
	}

	public void computeVectorField(double t) {
		xBuffer.position(0);
		xBuffer.get(x);
		ode.computeVectorField(t, x, xDot);
		xDotBuffer.position(0);
		xDotBuffer.put(xDot);
	}

	@Override
	public void computeVectorField(double t, double[] x, double[] xDot) {
		ode.computeVectorField(t, x, xDot);
	}

}
