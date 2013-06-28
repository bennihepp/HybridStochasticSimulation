package ch.ethz.khammash.ode;
public interface EventObserver {

	enum EventAction {
		STOP, CONTINUE,
	}

	void initialize(double t0, double[] x0, double t1);

	EventAction report(int eventIndex, double t, double[] x);

}
