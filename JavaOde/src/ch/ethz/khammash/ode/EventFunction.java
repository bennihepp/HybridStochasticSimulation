package ch.ethz.khammash.ode;

public interface EventFunction {

    int getNumberOfEventValues();

    void computeEventValues(double t, double[] x, double[] values);

}
