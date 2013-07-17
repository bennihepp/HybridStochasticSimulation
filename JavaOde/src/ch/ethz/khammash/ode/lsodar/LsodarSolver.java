package ch.ethz.khammash.ode.lsodar;

import ch.ethz.khammash.ode.EventFunction;
import ch.ethz.khammash.ode.EventObserver;
import ch.ethz.khammash.ode.FiniteTimepointProvider;
import ch.ethz.khammash.ode.Ode;
import ch.ethz.khammash.ode.Solver;
import ch.ethz.khammash.ode.StateObserver;
import ch.ethz.khammash.ode.TimepointProvider;

public class LsodarSolver implements Solver {

	private static LsodarSolver instance;

    static {
        // Load the shared library liblsodarjni
        System.loadLibrary("lsodarjni");
    }

    private boolean initialized = false;
    private long jni_pointer;
    private Ode ode;
    private EventFunction ef;
    private StateObserver stateObserver;
    private EventObserver eventObserver;
    private int eventIndex;
	private TimepointProvider timepointProvider;
	private double[] x;
//    private boolean[] eventFlags;
    private double relTol;
    private double absTol;

    public static synchronized LsodarSolver getInstance() {
    	if (instance == null)
            instance = new LsodarSolver();
    	return instance;
    }

    private LsodarSolver() {
    	relTol = 1e-6;
    	absTol = 1e-3;
    };

    private native long jni_initialize(Ode ode, EventFunction ef, double relTol, double absTol);

    public void setRelativeTolerance(double relTol) {
//    	checkArgument(relTol >= 0);
    	this.relTol = relTol;
    }

    public void setAbsoluteTolerance(double absTol) {
//    	checkArgument(absTol >= 0);
    	this.absTol = absTol;
    }

    public void initialize(Ode ode, EventFunction ef, StateObserver stateObserver, EventObserver eventObserver) throws JniInitializationException {
    	if (initialized)
    		dispose();
    	this.ode = ode;
    	this.ef = ef;
    	this.stateObserver = stateObserver;
    	this.eventObserver = eventObserver;
    	eventIndex = -1;
//    	eventFlags = new boolean[ode.getNumberOfEventValues()];
        jni_pointer = jni_initialize(ode, ef, relTol, absTol);
        if (jni_pointer == 0)
        	throw new JniInitializationException("Failed to initialize native LsodarDirectSolver");
        initialized = true;
    }

    private native void jni_dispose(long jni_pointer);

    public void dispose() {
    	if (initialized) {
    		jni_dispose(jni_pointer);
    		initialized = false;
    	}
    }

    private native double jni_integrate(long jni_pointer, double t0, double[] x0, double t1);

	@Override
	public double integrate(double t0, double[] x0, double t1) throws NotImplementedException {
		throw new NotImplementedException("Not implemented in LsodarSolver");
	}

	@Override
	public void prepare(TimepointProvider timepointProvider, double[] x) {
		this.timepointProvider = timepointProvider;
		this.x = x;
	}

    public double integrate() throws LsodarIntegrationException, NotYetInitializedException {
    	if (initialized) {
    		double t = timepointProvider.getCurrentTimepoint();
    		double tNext = timepointProvider.getNextTimepoint(t);
			stateObserver.report(t, x);
//    		double t0 = timepointProvider.getCurrentTimepoint();
//    		double t1 = timepointProvider.getLastTimepoint();
//    		stateObserver.initialize(t0, x, t1);
//    		eventObserver.initialize(t0, x, t1);
//    		double t = timepointProvider.getCurrentTimepoint();
//    		while (timepointProvider.hasNextTimepoint()) {
//    			double tNext = timepointProvider.getNextTimepoint();
    		while (true) {
//        		System.out.println("t="+t+", tNext="+tNext);
				t = jni_integrate(jni_pointer, t, x, tNext);
//        		System.out.println("prop=" + x0[x0.length-1]);
        		if (t < tNext)
            		if (eventOccured()) {
            			eventObserver.report(getEventOccuredIndex(), t, x);
            			resetEventOccuredFlags();
            			return t;
            		} else
            			throw new LsodarIntegrationException("Integration of ODE failed");
        		if (t == tNext)
        			stateObserver.report(t, x);
        		if (timepointProvider.hasNextTimepoint(t))
        			tNext = timepointProvider.getNextTimepoint(t);
        		else
        			break;
    		}
    		return t;
    	} else
    		throw new NotYetInitializedException("Solver has not yet been initialized");
    }

	@Override
	public double integrate(TimepointProvider timepointProvider, double[] x) throws LsodarIntegrationException {
		prepare(timepointProvider, x);
		return integrate();
	}

    private void reportEvent(int eventIndex, double t, double[] x) {
//    	eventFlags[eventIndex] = true;
    	if (this.eventIndex == -1)
    		this.eventIndex = eventIndex;
    }

    final private boolean eventOccured() {
    	return eventIndex != -1;
    }

    final private int getEventOccuredIndex() {
    	return eventIndex;
    }

//    final private boolean[] getEventOccuredFlags() {
//    	return eventFlags;
//    }

    final private void resetEventOccuredFlags() {
    	eventIndex = -1;
//    	Arrays.fill(eventFlags, false);
    }

    @Override
    protected void finalize() throws Throwable {
    	dispose();
    	super.finalize();
    }

    //private static native Object callbackTest2(double arr[]);

    //private static native int callbackTest(double arr[]);

    //public static double callback(double v) {
    //    return v + 1;
    //}

    public static void main(String args[]) {
    	LsodarSolver q = new LsodarSolver();
    	MyOde ode = new MyOde();
    	MyEf ef = new MyEf();
    	MyObserver observer = new MyObserver();
    	MyEventObserver eventObserver = new MyEventObserver();
    	double[] tSeries = { 0, 0.5, 1.0 };
    	q.initialize(ode, ef, observer, eventObserver);
    	FiniteTimepointProvider timepointProvider = new FiniteTimepointProvider(tSeries);
    	double[] x0 = { 0 };
    	q.integrate(timepointProvider, x0);
        q.dispose();
    };

    public static class MyOde implements Ode {

		@Override
		public int getDimensionOfVectorField() {
			return 1;
		}

		@Override
		public void computeVectorField(double t, double[] x, double[] xDot) {
//			System.out.println("computeVectorField");
			xDot[0] = x[0];
			xDot[0] = x[0] + 2;
		}

    }

    public static class MyEf implements EventFunction {

		@Override
		public int getNumberOfEventValues() {
			return 1;
		}

		@Override
		public void computeEventValues(double t, double[] x, double[] values) {
			values[0] = 2.0 - x[0];
		}
    	
    }

    public static class MyObserver implements StateObserver {

    	@Override
    	public void report(double t, double[] x) {
    		System.out.println("Report x["+t+"] = "+x[0]);
    	}

		@Override
		public void initialize(double t0, double[] x0, double t1) {
			report(t0, x0);
		}

    }

    public static class MyEventObserver implements EventObserver {

		@Override
		public void initialize(double t0, double[] x0, double t1) {
		}

		@Override
		public EventAction report(int eventIndex, double t, double[] x) {
	        System.out.println("Event "+eventIndex+" occured at t="+t);
	        return EventAction.STOP;
		}

    }

}
