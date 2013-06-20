package ch.ethz.khammash.nativeode;

public interface EventFunction {

    int getNumberOfEventValues();

    void computeEventValues(double t, double[] x, double[] values);

}
