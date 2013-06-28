package ch.ethz.khammash.ode;
public interface StateObserver {

	void initialize(double t0, double[] x0, double t1);

	void report(double t, double[] x);

}
