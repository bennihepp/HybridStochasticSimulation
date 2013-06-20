package ch.ethz.khammash.nativeode;

public interface Solver {

    void initialize(Ode ode, EventFunction ef, StateObserver stateObserver, EventObserver eventObserver);

    void dispose();

    public double integrate(TimepointIterator timepointIterator, double[] x0);

}
