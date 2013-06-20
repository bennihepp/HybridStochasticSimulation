package ch.ethz.khammash.nativeode;
public interface StateObserver {

	void initialize(double t0, double[] x0, double t1);

	void report(double t, double[] x);

}
