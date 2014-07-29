package ch.ethz.bhepp.ode.cvode;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;

import ch.ethz.bhepp.ode.AdaptiveStepSolver;
import ch.ethz.bhepp.ode.BufferEventFunctionAdapter;
import ch.ethz.bhepp.ode.BufferOdeAdapter;
import ch.ethz.bhepp.ode.DummyEventFunction;
import ch.ethz.bhepp.ode.DummyEventObserver;
import ch.ethz.bhepp.ode.EventFunction;
import ch.ethz.bhepp.ode.EventObserver;
import ch.ethz.bhepp.ode.FiniteTimepointProvider;
import ch.ethz.bhepp.ode.Ode;
import ch.ethz.bhepp.ode.Solver;
import ch.ethz.bhepp.ode.StateObserver;
import ch.ethz.bhepp.ode.TimepointProvider;
import ch.ethz.bhepp.ode.EventObserver.EventAction;

public class CVodeSolver implements Solver, AdaptiveStepSolver {

	private static volatile boolean cvodeJNILibraryLoaded = false;

    private static void ensureCVodeJNILibraryLoaded() {
    	if (!cvodeJNILibraryLoaded) {
    		synchronized (CVodeSolver.class) {
    			if (!cvodeJNILibraryLoaded) {
    				loadCVodeJNILibrary();
    				cvodeJNILibraryLoaded = true;
    			}
    		}
    	}
	}

    private static void loadCVodeJNILibrary() {
        // Load the shared library libcvodejni
        System.loadLibrary("cvodejni");
    }

	public static class JniException extends RuntimeException {
		private static final long serialVersionUID = 4397470726754969807L;
		private int errorCode;

		public JniException(String message) {
			super(message);
			errorCode = 0;
		}

		public JniException(String message, int errorCode) {
			super(message + " (error code=" + errorCode + ")");
			this.errorCode = errorCode;
		}

		public int getErrorCode() {
			return errorCode;
		}
	}

    public static final int MULTISTEPTYPE_ADAMS = 1;
    public static final int MULTISTEPTYPE_BDF = 2;
    public static final int ITERATIONTYPE_FUNCTIONAL = 1;
    public static final int ITERATIONTYPE_NEWTON = 2;

    private boolean initialized = false;
    private long jni_pointer;
    private BufferOdeAdapter ode;
    private BufferEventFunctionAdapter ef;
    private StateObserver stateObserver;
    private EventObserver eventObserver;
	private TimepointProvider timepointProvider;
	private double[] x;
	private double t1;
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
	//private IntBuffer eventIndexBuffer;
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
    	ensureCVodeJNILibraryLoaded();
    	this.relTolerance = relTolerance;
    	this.absTolerance = absTolerance;
    };

	private native long jni_initialize(BufferOdeAdapter ode, BufferEventFunctionAdapter ef,
			DoubleBuffer xBuffer, DoubleBuffer xTmpBuffer, DoubleBuffer xDotBuffer,
			DoubleBuffer gBuffer,
			double relTol, double absTol, int multistepType, int iterationType, int maxNumOfSteps, double minStep, double maxStep);

    private native void jni_dispose(long jni_pointer);

    private native void jni_resetEventIndex(long jni_pointer);

    private native int jni_getEventIndex(long jni_pointer);

    private native void jni_reinitializeOneStep(long jni_pointer, double t);

    private native double jni_integrateOneStep(long jni_pointer, double t1);

    private native void jni_computeInterpolatedSolution(long jni_pointer, double t1);

    private native double jni_getCurrentStep(long jni_pointer);

    private native double jni_getLastStep(long jni_pointer);

    private native void jni_setInitStep(long jni_pointer, double step);

    private native void jni_setMinStep(long jni_pointer, double minStep);

    private native void jni_setMaxStep(long jni_pointer, double maxStep);

    private native void jni_setStopTime(long jni_pointer, double stopTime);

    private native void jni_reinitialize(long jni_pointer, double t);

    private native void jni_quick_reinitialize(long jni_pointer, double t);

    private native double jni_integrate(long jni_pointer, double t1);

	public void setMultistepType(int multistepType) throws IllegalArgumentException {
    	if (multistepType != MULTISTEPTYPE_ADAMS && multistepType != MULTISTEPTYPE_BDF)
    		throw new IllegalArgumentException("Multistep type must be either MULTISTEPTYPE_ADAMS or MULTISTEPTYPE_BDF");
    	this.multistepType = multistepType;
    }

    public void setIterationType(int iterationType) throws IllegalArgumentException {
    	if (iterationType != ITERATIONTYPE_FUNCTIONAL && iterationType != ITERATIONTYPE_NEWTON)
    		throw new IllegalArgumentException("Iteration type must be either ITERATIONTYPE_FUNCTIONAL or ITERATIONTYPE_NEWTON");
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

	@Override
	public void initialize(Ode ode) throws InitializationException {
		initialize(ode, null, null, null);
	}

	@Override
    public void initialize(Ode ode, EventFunction ef, StateObserver stateObserver, EventObserver eventObserver) throws JniException {
    	if (ef == null) {
    		ef = new DummyEventFunction();
    		eventObserver = new DummyEventObserver();
    	}
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
    	//eventIndexBuffer = ibuffer.slice();

    	this.ode = new BufferOdeAdapter(ode, xTmpBuffer, xDotBuffer);
    	this.ef = new BufferEventFunctionAdapter(ef, ode.getDimensionOfVectorField(), xTmpBuffer, gBuffer);
    	try {
    		jni_pointer = jni_initialize(this.ode, this.ef, xBuffer, xTmpBuffer, xDotBuffer, gBuffer,
    				relTolerance, absTolerance, multistepType, iterationType, maxNumOfSteps, minStep, maxStep);
    	} catch (JniException e) {
    		throw new InitializationException("Failed to initialize native CVodeSolver", e);
    	}
        if (jni_pointer == 0)
        	throw new JniException("Failed to initialize native CVodeSolver");
        initialized = true;
		prepared = false;
    }

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

	@Override
	public double integrate(double t0, double[] x0, double t1) throws UnsupportedOperationException {
		throw new UnsupportedOperationException("Not implemented in CVodeSolver");
	}

	@Override
	public void prepareStep(double t0, double[] x0, double t1) throws IllegalStateException {
    	if (!initialized)
    		throw new IllegalStateException("Solver has not yet been initialized");

		this.x = x0;
		xBuffer.position(0);
		xBuffer.put(x);
		jni_reinitializeOneStep(jni_pointer, t0);
		resetEventOccuredFlags();
		this.t1 = t1;
	}

	@Override
	public void prepareStep(double t0, double[] x0, double t1, double dt0) throws IllegalStateException {
		prepareStep(t0, x0, t1);
		setCurrentStepSize(dt0);
	}

	@Override
	public double integrateStep() {
		double tRet = jni_integrateOneStep(jni_pointer, t1);
		xBuffer.position(0);
		xBuffer.get(x);

//		if (eventOccured()) {
//			eventObserver.report(getEventOccuredIndex(), tRet, x);
//			resetEventOccuredFlags();
//		}

		return tRet;
	}

	@Override
	public double integrateStep(double[] xOut) {
		double tRet = jni_integrateOneStep(jni_pointer, t1);
		xBuffer.position(0);
		xBuffer.get(xOut);

//		if (eventOccured()) {
//			eventObserver.report(getEventOccuredIndex(), tRet, x);
//			resetEventOccuredFlags();
//		}

		return tRet;
	}

	@Override
	public void computeInterpolatedSolution(double t, double[] xOut) throws IllegalStateException {
		jni_computeInterpolatedSolution(jni_pointer, t);
		xTmpBuffer.position(0);
		xTmpBuffer.get(xOut);
	}

//	public double integrateOneStep(double tNext) throws CVodeIntegrationException, NotYetInitializedException {
////    	if (!initialized)
////    		throw new NotYetInitializedException("Solver has not yet been initialized");
//
//		double tRet = jni_integrateOneStep(jni_pointer, tNext);
//		xBuffer.position(0);
//		xBuffer.get(x);
//
//		if (eventOccured()) {
//			eventObserver.report(getEventOccuredIndex(), tRet, x);
//			resetEventOccuredFlags();
//		}
//
//		return tRet;
//	}

	public void prepare(double t, double[] x) {
		this.x = x;
		xBuffer.position(0);
		xBuffer.put(x);
		jni_reinitialize(jni_pointer, t);
		resetEventOccuredFlags();
	}

	public double integrate(double tNext) throws IntegrationException, IllegalStateException {
//    	if (!initialized)
//    		throw new NotYetInitializedException("Solver has not yet been initialized");
		double tRet = jni_integrate(jni_pointer, tNext);
		xBuffer.position(0);
		xBuffer.get(x);

		if (tRet < tNext)
    		if (eventOccured()) {
    			EventAction ea = eventObserver.report(getEventOccuredIndex(), tRet, x);
    			resetEventOccuredFlags();
    			if (ea == EventAction.STOP)
        			return tRet;
    		} else
    			throw new IntegrationException("Integration of ODE failed: tRet < tNext");
		if (tRet > tNext)
			throw new RuntimeException("Integration of ODE failed: tRet > tNext");

		return tRet;
	}

	@Override
	public void prepare(TimepointProvider timepointProvider, double[] x) {
		this.timepointProvider = timepointProvider;
		this.x = x;
		xBuffer.position(0);
		xBuffer.put(x);
		double t = timepointProvider.getCurrentTimepoint();
		jni_reinitialize(jni_pointer, t);
		resetEventOccuredFlags();
	}

	@Override
	public double integrate() throws IntegrationException, IllegalStateException {
    	if (!initialized)
    		throw new IllegalStateException("Solver has not yet been initialized");

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
			t = jni_integrate(jni_pointer, tNext);
    		xBuffer.position(0);
    		xBuffer.get(x);

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
        			throw new IntegrationException("Integration of ODE failed: t < tNext");
    		if (t > tNext)
    			throw new RuntimeException("Integration of ODE failed: t > tNext");
			stateObserver.report(t, x);
		}
		return t;
	}

    public double integrate(TimepointProvider timepointProvider, double[] x) throws IntegrationException {
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
//    	return eventIndexBuffer.get(0) != -1;
    	return getEventOccuredIndex() != -1;
    }

	final private int getEventOccuredIndex() {
//    	return eventIndex;
//    	return eventIndexBuffer.get(0);
    	return jni_getEventIndex(jni_pointer);
    }

//    final private boolean[] getEventOccuredFlags() {
//    	return eventFlags;
//    }

    final private void resetEventOccuredFlags() {
//    	eventIndex = -1;
//    	eventIndexBuffer.put(0, -1);
//    	Arrays.fill(eventFlags, false);
    	jni_resetEventIndex(jni_pointer);
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

	@Override
	public double getCurrentStepSize() throws IllegalStateException, UnsupportedOperationException {
		return jni_getCurrentStep(jni_pointer);
	}

	@Override
	public void setCurrentStepSize(double stepSize) throws IllegalStateException, UnsupportedOperationException {
		jni_setInitStep(jni_pointer, stepSize);
	}

	@Override
	public double getLastStepSize() throws IllegalStateException, UnsupportedOperationException {
		return jni_getLastStep(jni_pointer);
	}

	@Override
	public void setMinimumStepSize(double minStepSize) throws IllegalStateException, UnsupportedOperationException {
		jni_setMinStep(jni_pointer, minStepSize);
	}

	@Override
	public void setMaximumStepSize(double maxStepSize) throws IllegalStateException, UnsupportedOperationException {
		jni_setMaxStep(jni_pointer, maxStepSize);
	}

    public void setStopTime(double stopTime) {
    	jni_setStopTime(jni_pointer, stopTime);
    }

}
