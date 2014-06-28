package ch.ethz.bhepp.ode;

public class DummyEventFunction implements EventFunction {

	@Override
	public int getNumberOfEventValues() {
		return 0;
	}

	@Override
	public void computeEventValues(double t, double[] x, double[] values) {
	}

}
