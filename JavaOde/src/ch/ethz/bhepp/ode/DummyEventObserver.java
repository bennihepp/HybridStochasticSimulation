package ch.ethz.bhepp.ode;

public class DummyEventObserver implements EventObserver {

	@Override
	public void initialize(double t0, double[] x0, double t1) {
	}

	@Override
	public EventAction report(int eventIndex, double t, double[] x) {
		return EventAction.CONTINUE;
	}

}
