package ch.ethz.khammash.ode.lsodar;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;

import ch.ethz.khammash.ode.BufferEventFunctionAdapter;
import ch.ethz.khammash.ode.BufferOdeAdapter;
import ch.ethz.khammash.ode.EventFunction;
import ch.ethz.khammash.ode.EventObserver;
import ch.ethz.khammash.ode.FiniteTimepointProvider;
import ch.ethz.khammash.ode.Ode;
import ch.ethz.khammash.ode.Solver;
import ch.ethz.khammash.ode.StateObserver;
import ch.ethz.khammash.ode.TimepointProvider;

public class LsodarDirectSolver implements Solver {

	private static LsodarDirectSolver instance;

    static {
        // Load the shared library liblsodarjni
        System.loadLibrary("lsodarjni");
    }

    private boolean initialized = false;
    private long jni_pointer;
    private BufferOdeAdapter lsodarOde;
    private BufferEventFunctionAdapter lsodarEf;
    private StateObserver stateObserver;
    private EventObserver eventObserver;
	private TimepointProvider timepointProvider;
	private double[] x;
//    private int eventIndex;
//    private boolean[] eventFlags;
    private double relTol;
    private double absTol;
	private DoubleBuffer xBuffer;
	private DoubleBuffer xTmpBuffer;
	private DoubleBuffer xDotBuffer;
	private DoubleBuffer gBuffer;
	private IntBuffer eventIndexBuffer;

    public static synchronized LsodarDirectSolver getInstance() {
    	if (instance == null)
            instance = new LsodarDirectSolver();
    	return instance;
    }

    private LsodarDirectSolver() {
    	relTol = 1e-6;
    	absTol = 1e-3;
    };

//    public double[] getxArray() {
//    	return x;
//    }
//
//    public double[] getxTmpArray() {
//    	return xTmp;
//    }
//
//    public double[] getxDotArray() {
//    	return xDot;
//    }
//
//    public double[] getgArray() {
//    	return g;
//    }

    public void setRelativeTolerance(double relTol) {
//    	checkArgument(relTol >= 0);
    	this.relTol = relTol;
    }

    public void setAbsoluteTolerance(double absTol) {
//    	checkArgument(absTol >= 0);
    	this.absTol = absTol;
    }

	private native long jni_initialize(BufferOdeAdapter ode, BufferEventFunctionAdapter ef,
			DoubleBuffer xBuffer, DoubleBuffer xTmpBuffer, DoubleBuffer xDotBuffer,
			DoubleBuffer gBuffer, IntBuffer eventIndexBuffer,
			double relTol, double absTol);

    public void initialize(Ode ode, EventFunction ef, StateObserver stateObserver, EventObserver eventObserver) throws JniInitializationException {
    	if (initialized)
    		dispose();
    	this.stateObserver = stateObserver;
    	this.eventObserver = eventObserver;
//    	eventIndex = -1;
//    	eventFlags = new boolean[ode.getNumberOfEventValues()];

    	int bufferSize = Double.SIZE / 8 * (3 * ode.getDimensionOfVectorField() + ef.getNumberOfEventValues())
    					+ Integer.SIZE / 8;
    	ByteBuffer buffer = ByteBuffer.allocateDirect(bufferSize);
    	buffer.order(ByteOrder.nativeOrder());
    	DoubleBuffer dbuffer = buffer.asDoubleBuffer();
    	IntBuffer ibuffer = buffer.asIntBuffer();
    	int offset = 0;
    	dbuffer.position(offset);
    	xBuffer = dbuffer.slice();
    	offset += ode.getDimensionOfVectorField();
    	dbuffer.position(offset);
    	xTmpBuffer = dbuffer.slice();
    	offset += ode.getDimensionOfVectorField();
    	dbuffer.position(offset);
    	xDotBuffer = dbuffer.slice();
    	offset += ode.getDimensionOfVectorField();
    	dbuffer.position(offset);
    	gBuffer = dbuffer.slice();
    	offset += ef.getNumberOfEventValues();
    	ibuffer.position((Double.SIZE / Integer.SIZE) / 8 * offset);
    	eventIndexBuffer = ibuffer.slice();

//    	xDot = xDotBuffer.asDoubleBuffer().array();
//    	x = xBuffer.asDoubleBuffer().array();
//    	xTmp = xTmpBuffer.asDoubleBuffer().array();
//    	g = gBuffer.asDoubleBuffer().array();
//    	eventIndex = eventIndexBuffer.asIntBuffer().array();

//    	xBuffer = ByteBuffer.allocateDirect(Double.SIZE / 8 * ode.getDimensionOfVectorField());
//    	x = xBuffer.asDoubleBuffer().array();
//    	xTmpBuffer = ByteBuffer.allocateDirect(Double.SIZE / 8 * ode.getDimensionOfVectorField());
//    	xTmp = xTmpBuffer.asDoubleBuffer().array();
//    	xDotBuffer = ByteBuffer.allocateDirect(Double.SIZE / 8 * ode.getDimensionOfVectorField());
//    	xDot = xDotBuffer.asDoubleBuffer().array();
//    	gBuffer = ByteBuffer.allocateDirect(Double.SIZE / 8 * ef.getNumberOfEventValues());
//    	g = gBuffer.asDoubleBuffer().array();

    	this.lsodarOde = new BufferOdeAdapter(ode, xTmpBuffer, xDotBuffer);
    	this.lsodarEf = new BufferEventFunctionAdapter(ef, ode.getDimensionOfVectorField(), xTmpBuffer, gBuffer);
        jni_pointer = jni_initialize(lsodarOde, lsodarEf, xBuffer, xTmpBuffer, xDotBuffer, gBuffer, eventIndexBuffer, relTol, absTol);
//        jni_pointer = jni_initialize(lsodarOde, lsodarEf, xBuffer, xTmpBuffer, xDotBuffer, gBuffer, eventIndexBuffer, relTol, absTol);
        if (jni_pointer == 0)
        	throw new JniInitializationException("Failed to initialize native LsodarDirectSolver");
        initialized = true;
    }

    private native void jni_dispose(long jni_pointer);

    public void dispose() {
    	if (initialized) {
//    		x = null;
    		xBuffer = null;
//    		xTmp = null;
    		xTmpBuffer = null;
//    		xDot = null;
    		xDotBuffer = null;
//    		g = null;
    		gBuffer = null;
    		jni_dispose(jni_pointer);
    		initialized = false;
    	}
    }

	private native void jni_reinitialize(long jni_pointer, double t);

    private native double jni_integrate(long jni_pointer, double t0, double t1);

	@Override
	public void prepare(TimepointProvider timepointProvider, double[] x) {
		this.timepointProvider = timepointProvider;
		this.x = x;
	}

	@Override
	public double integrate(double t0, double[] x0, double t1) {
		throw new NotImplementedException("Not implemented in LsodarDirectSolver");
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
        		xBuffer.put(x);
        		xBuffer.position(0);
				t = jni_integrate(jni_pointer, t, tNext);
	    		xBuffer.get(x);
        		xBuffer.position(0);
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

//    private void reportEvent(int eventIndex, double t, double[] x) {
////    	eventFlags[eventIndex] = true;
//    	if (this.eventIndex == -1)
//    		this.eventIndex = eventIndex;
//    }

    final private boolean eventOccured() {
//    	return eventIndex != -1;
    	return eventIndexBuffer.get(0) != -1;
    }

    final private int getEventOccuredIndex() {
//    	return eventIndex;
    	return eventIndexBuffer.get(0);
    }

//    final private boolean[] getEventOccuredFlags() {
//    	return eventFlags;
//    }

    final private void resetEventOccuredFlags() {
//    	eventIndex = -1;
    	eventIndexBuffer.put(0, -1);
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
    	LsodarDirectSolver q = new LsodarDirectSolver();
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
