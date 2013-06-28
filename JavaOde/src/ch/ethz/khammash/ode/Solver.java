package ch.ethz.khammash.ode;

public interface Solver {

    void initialize(Ode ode, EventFunction ef, StateObserver stateObserver, EventObserver eventObserver);

    void dispose();

    public void prepare(TimepointProvider timepointProvider, double[] x);

    // User is responsible to call prepare first
    public double integrate();

    public double integrate(double t0, double[] x0, double t1);

    public double integrate(TimepointProvider timepointProvider, double[] x0);

}
