package test;

import ch.ethz.bhepp.ode.Ode;


public class SimpleOde implements Ode {

	@Override
	public int getDimensionOfVectorField() {
		return 1;
	}

	@Override
	public void computeVectorField(double t, double[] x, double[] xDot) {
		xDot[0] = x[0];
		xDot[0] = x[0] + 2;
	}

}
