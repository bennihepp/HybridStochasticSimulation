package ch.ethz.khammash.ode.cvode;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;

import ch.ethz.khammash.ode.DirectBufferEventFunctionAdapter;
import ch.ethz.khammash.ode.DirectBufferOdeAdapter;
import ch.ethz.khammash.ode.EventFunction;
import ch.ethz.khammash.ode.EventObserver;
import ch.ethz.khammash.ode.EventObserver.EventAction;
import ch.ethz.khammash.ode.Exception;
import ch.ethz.khammash.ode.FiniteTimepointProvider;
import ch.ethz.khammash.ode.Ode;
import ch.ethz.khammash.ode.Solver;
import ch.ethz.khammash.ode.StateObserver;
import ch.ethz.khammash.ode.TimepointProvider;

public class CVodeSolver implements Solver {

	public static class InvalidMultistepTypeException extends Exception {
		private static final long serialVersionUID = 691615125296233066L;

		public InvalidMultistepTypeException(String message) {
			super(message);
		}
	}

	public static class InvalidIterationTypeException extends Exception {
		private static final long serialVersionUID = -742731064807994470L;

		public InvalidIterationTypeException(String message) {
			super(message);
		}
	}

	public static class JniInitializationException extends InitializationException {
		private static final long serialVersionUID = -742731064807994470L;

		public JniInitializationException(String message) {
			super(message);
		}
	}

	static class CVodeIntegrationException extends IntegrationException {
		private static final long serialVersionUID = 8509338706972523801L;

		public CVodeIntegrationException(String message) {
			super(message);
		}

		public CVodeIntegrationException(String message, Throwable cause) {
			super(message, cause);
		}
	}

    static {
        // Load the shared library libcvodejni
        System.loadLibrary("cvodejni");
    }

    public static int MULTISTEPTYPE_ADAMS = 1;
    public static int MULTISTEPTYPE_BDF = 2;
    public static int ITERATIONTYPE_FUNCTIONAL = 1;
    public static int ITERATIONTYPE_NEWTON = 2;

    private boolean initialized = false;
    private long jni_pointer;
    private DirectBufferOdeAdapter ode;
    private DirectBufferEventFunctionAdapter ef;
    private StateObserver stateObserver;
    private EventObserver eventObserver;
	private TimepointProvider timepointProvider;
	private double[] x;
	protected double[] eventValues;
	protected int[] eventSigns;
//    private int eventIndex;
//    private boolean[] eventFlags;
    private double relTolerance;
    private double absTolerance;
	private DoubleBuffer xBuffer;
	private DoubleBuffer xTmpBuffer;
	private DoubleBuffer xDotBuffer;
	private DoubleBuffer gBuffer;
	private IntBuffer eventIndexBuffer;
	private boolean prepared;
	private int multistepType = MULTISTEPTYPE_BDF;
	private int iterationType = ITERATIONTYPE_NEWTON;
	private int maxNumOfSteps = 0;
	private double minStep = Double.NaN;
	private double maxStep = Double.NaN;

    public CVodeSolver() {
    	this(1e-6, 1e-3);
    };

    public CVodeSolver(double relTolerance, double absTolerance) {
    	this.relTolerance = relTolerance;
    	this.absTolerance = absTolerance;
    };

    public void setMultistepType(int multistepType) throws InvalidMultistepTypeException {
    	if (multistepType != MULTISTEPTYPE_ADAMS && multistepType != MULTISTEPTYPE_BDF)
    		throw new InvalidMultistepTypeException("Multistep type must be either MULTISTEPTYPE_ADAMS or MULTISTEPTYPE_BDF");
    	this.multistepType = multistepType;
    }

    public void setIterationType(int iterationType) throws InvalidIterationTypeException {
    	if (iterationType != ITERATIONTYPE_FUNCTIONAL && iterationType != ITERATIONTYPE_NEWTON)
    		throw new InvalidIterationTypeException("Iteration type must be either ITERATIONTYPE_FUNCTIONAL or ITERATIONTYPE_NEWTON");
    	this.iterationType = iterationType;
    }

    public void setRelativeTolerance(double relTolerance) {
//    	checkArgument(relTol >= 0);
    	this.relTolerance = relTolerance;
    }

    public void setAbsoluteTolerance(double absTolerance) {
//    	checkArgument(absTol >= 0);
    	this.absTolerance = absTolerance;
    }

    public void setMaxNumOfSteps(int maxNumOfSteps) {
    	this.maxNumOfSteps = maxNumOfSteps;
    }

    public void setMinStep(double minStep) {
    	this.minStep = minStep;
    }

    public void setMaxStep(double maxStep) {
    	this.maxStep = maxStep;
    }

	private native long jni_initialize(DirectBufferOdeAdapter ode, DirectBufferEventFunctionAdapter ef,
			DoubleBuffer xBuffer, DoubleBuffer xTmpBuffer, DoubleBuffer xDotBuffer,
			DoubleBuffer gBuffer, IntBuffer eventIndexBuffer,
			double relTol, double absTol, int multistepType, int iterationType, int maxNumOfSteps, double minStep, double maxStep);

	@Override
    public void initialize(Ode ode, EventFunction ef, StateObserver stateObserver, EventObserver eventObserver) throws JniInitializationException {
    	if (initialized)
    		dispose();
    	this.stateObserver = stateObserver;
    	this.eventObserver = eventObserver;
		this.eventValues = new double[ef.getNumberOfEventValues()];
		this.eventSigns = new int[ef.getNumberOfEventValues()];
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

    	this.ode = new DirectBufferOdeAdapter(ode, xTmpBuffer, xDotBuffer);
    	this.ef = new DirectBufferEventFunctionAdapter(ef, ode.getDimensionOfVectorField(), xTmpBuffer, gBuffer);
        jni_pointer = jni_initialize(this.ode, this.ef, xBuffer, xTmpBuffer, xDotBuffer, gBuffer, eventIndexBuffer,
        		relTolerance, absTolerance, multistepType, iterationType, maxNumOfSteps, minStep, maxStep);
        if (jni_pointer == 0)
        	throw new JniInitializationException("Failed to initialize native CVodeSolver");
        initialized = true;
		prepared = false;
    }

    private native void jni_dispose(long jni_pointer);

    public void dispose() {
    	if (initialized) {
    		xBuffer = null;
    		xTmpBuffer = null;
    		xDotBuffer = null;
    		gBuffer = null;
    		jni_dispose(jni_pointer);
    		initialized = false;
    	}
    }

    private native void jni_reinitialize(long jni_pointer, double t);

    private native void jni_quick_reinitialize(long jni_pointer, double t);


    private native double jni_integrate(long jni_pointer, double t1);

	@Override
	public double integrate(double t0, double[] x0, double t1) throws NotImplementedException {
		throw new NotImplementedException("Not implemented in CVodeSolver");
	}

	@Override
	public void prepare(TimepointProvider timepointProvider, double[] x) {
		this.timepointProvider = timepointProvider;
		this.x = x;
		xBuffer.position(0);
		xBuffer.put(x);
		double t = timepointProvider.getCurrentTimepoint();
//		if (prepared)
//			jni_quick_reinitialize(jni_pointer, t);
//		else {
//			jni_reinitialize(jni_pointer, t);
//			prepared = true;
//		}
		jni_reinitialize(jni_pointer, t);
	}

	@Override
	public double integrate() throws CVodeIntegrationException, NotYetInitializedException {
    	if (initialized) {
    		double t = timepointProvider.getCurrentTimepoint();
			stateObserver.report(t, x);

//			// Crude root-finding
//			ef.computeEventValues(t, x, eventValues);
//			for (int i=0; i < ef.getNumberOfEventValues(); i++) {
//				if (eventValues[i] == 0)
//					throw new Exception("Event function has zeros at start of integration");
//				eventSigns[i] = (eventValues[i] > 0) ? 1 : -1;
//			}

//    		double t0 = timepointProvider.getCurrentTimepoint();
//    		double t1 = timepointProvider.getLastTimepoint();
//    		stateObserver.initialize(t0, x, t1);
//    		eventObserver.initialize(t0, x, t1);
//    		double t = timepointProvider.getCurrentTimepoint();
//    		while (timepointProvider.hasNextTimepoint()) {
//    			double tNext = timepointProvider.getNextTimepoint();
			while (timepointProvider.hasNextTimepoint(t)) {
        		double tNext = timepointProvider.getNextTimepoint(t);
//        		System.out.println("t="+t+", tNext="+tNext);
				t = jni_integrate(jni_pointer, tNext);
        		xBuffer.position(0);
	    		xBuffer.get(x);
//        		System.out.println("prop=" + x0[x0.length-1]);

//				// Crude root-finding
//				ef.computeEventValues(t, x, eventValues);
//				for (int k=0; k < ef.getNumberOfEventValues(); k++)
//					if (eventSigns[k] * eventValues[k] < 0) {
//						// Event occurred during last step
//						if (eventObserver.report(k, t, x) == EventAction.STOP)
//							return t;
//	        		}

        		if (t < tNext)
            		if (eventOccured()) {
            			EventAction ea = eventObserver.report(getEventOccuredIndex(), t, x);
            			resetEventOccuredFlags();
            			if (ea == EventAction.STOP)
	            			return t;
            		} else
            			throw new CVodeIntegrationException("Integration of ODE failed: t < tNext");
        		if (t > tNext)
        			throw new RuntimeException("Integration of ODE failed: t > tNext");
    			stateObserver.report(t, x);
    		}
    		return t;
    	} else
    		throw new NotYetInitializedException("Solver has not yet been initialized");
	}

    public double integrate(TimepointProvider timepointProvider, double[] x) throws CVodeIntegrationException {
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

    public static void main(String args[]) {
    	CVodeSolver q = new CVodeSolver();
    	MyOde ode = new MyOde();
    	MyEf ef = new MyEf();
    	MyObserver observer = new MyObserver();
    	MyEventObserver eventObserver = new MyEventObserver();
    	double[] tSeries = { 0, 0.5, 1.0 };
    	q.initialize(ode, ef, observer, eventObserver);
    	FiniteTimepointProvider timepointProvider = new FiniteTimepointProvider(tSeries);
    	double[] x0 = { 0.2 };
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
