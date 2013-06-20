package ch.ethz.khammash.nativeode;
public interface EventObserver {

	void initialize(double t0, double[] x0, double t1);

	void report(int eventIndex, double t, double[] x);

}
