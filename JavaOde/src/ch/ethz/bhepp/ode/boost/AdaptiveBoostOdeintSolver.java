package ch.ethz.bhepp.ode.boost;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;

import org.apache.commons.math3.util.FastMath;

import ch.ethz.bhepp.ode.AdaptiveStepSolver;
import ch.ethz.bhepp.ode.BufferOdeAdapter;
import ch.ethz.bhepp.ode.Ode;
import ch.ethz.bhepp.ode.Solver.InitializationException;

public class AdaptiveBoostOdeintSolver implements AdaptiveStepSolver {

	private static volatile boolean boostJNILibraryLoaded = false;

    private static void ensureBoostJNILibraryLoaded() {
    	if (!boostJNILibraryLoaded) {
    		synchronized (FixedBoostOdeintSolver.class) {
    			if (!boostJNILibraryLoaded) {
    				loadBoostJNILibrary();
    				boostJNILibraryLoaded = true;
    			}
    		}
    	}
	}

    private static void loadBoostJNILibrary() {
        // Load the shared library libboostjni
        System.loadLibrary("boostjni");
    }

	private BufferOdeAdapter ode;
	private StepperType stepperType;
	private boolean initialized;
	private boolean prepared;
	private DoubleBuffer xBuffer;
	private DoubleBuffer xTmpBuffer;
	private DoubleBuffer xDotBuffer;
	private long jni_pointer;
	private double relTolerance;
	private double absTolerance;
	private double initialStepSize;
	private double t;
	private double t1;


	public AdaptiveBoostOdeintSolver(double initialStepSize, StepperType type) {
		ensureBoostJNILibraryLoaded();
		initialized = false;
		prepared = false;
		this.initialStepSize = initialStepSize;
		this.stepperType = type;
		relTolerance = 1e-3;
		absTolerance = 1e-3;
	}

    public void setRelativeTolerance(double relTolerance) {
//    	checkArgument(relTol >= 0);
    	this.relTolerance = relTolerance;
    }

    public void setAbsoluteTolerance(double absTolerance) {
//    	checkArgument(absTol >= 0);
    	this.absTolerance = absTolerance;
    }

	private native long jni_initialize(BufferOdeAdapter ode,
			DoubleBuffer xBuffer, DoubleBuffer xTmpBuffer, DoubleBuffer xDotBuffer,
			double initialStepSize, double relTol, double absTol,
			int stepperType);

	@Override
	public void initialize(Ode ode) throws InitializationException {
    	if (initialized)
    		dispose();

    	int bufferSize = Double.SIZE / 8 * (3 * ode.getDimensionOfVectorField())
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
    	ibuffer.position((Double.SIZE / Integer.SIZE) / 8 * offset);

    	this.ode = new BufferOdeAdapter(ode, xTmpBuffer, xDotBuffer);
    	try {
    		jni_pointer = jni_initialize(this.ode, xBuffer, xTmpBuffer, xDotBuffer,
    				initialStepSize, relTolerance, absTolerance, stepperType.getType());
    	} catch (JniException e) {
    		throw new InitializationException("Failed to initialize native CVodeSolver", e);
    	}
        if (jni_pointer == 0)
        	throw new JniException("Failed to initialize native CVodeSolver");
        initialized = true;
		prepared = false;
	}

    private native void jni_dispose(long jni_pointer);

	public void dispose() {
    	if (initialized) {
    		xBuffer = null;
    		xTmpBuffer = null;
    		xDotBuffer = null;
    		jni_dispose(jni_pointer);
    		initialized = false;
    	}
	}

    private native void jni_prepareStep(long jni_pointer, double t0, double t1);

	@Override
	public void prepareStep(double t0, double[] x0, double t1) throws IllegalStateException {
		prepared = true;
		this.t = t0;
		this.t1 = t1;
		xBuffer.position(0);
		xBuffer.put(x0);
		jni_prepareStep(jni_pointer, t0, t1);
	}

	private native void jni_setCurrentStepSize(long jni_pointer, double stepSize);

	@Override
	public void prepareStep(double t0, double[] x0, double t1, double dt0) throws IllegalStateException {
		jni_setCurrentStepSize(jni_pointer, dt0);
		prepareStep(t0, x0, t1);
	}

    private native double jni_getCurrentState(long jni_pointer);

	public double getCurrentState(double[] xOut) {
		double t = jni_getCurrentState(jni_pointer);
		xTmpBuffer.position(0);
		xTmpBuffer.get(xOut);
		return t;
	}

    private native double jni_integrateStep(long jni_pointer);

	@Override
	public double integrateStep() {
		assert(prepared);
		double tNext = jni_integrateStep(jni_pointer);
		return FastMath.min(tNext, t1);
	}

	@Override
	public double integrateStep(double[] xOut) {
		assert(prepared);
		double tNext = jni_integrateStep(jni_pointer);
		if (tNext > t1) {
			tNext = t1;
			computeInterpolatedSolution(tNext, xOut);
		} else {
			getCurrentState(xOut);
		}
		return tNext;
	}

	private native void jni_computeInterpolatedSolution(long jni_pointer, double t);

	@Override
	public void computeInterpolatedSolution(double t, double[] xOut) throws IllegalStateException {
		assert(prepared);
		jni_computeInterpolatedSolution(jni_pointer, t);
		xTmpBuffer.position(0);
		xTmpBuffer.get(xOut);
	}

	private native double jni_getCurrentStepSize(long jni_pointer);

	@Override
	public double getCurrentStepSize() throws IllegalStateException, UnsupportedOperationException {
		return jni_getCurrentStepSize(jni_pointer);
	}

	private native double jni_getLastStepSize(long jni_pointer);

	@Override
	public double getLastStepSize() throws IllegalStateException, UnsupportedOperationException {
		return jni_getLastStepSize(jni_pointer);
	}

	@Override
	public void setCurrentStepSize(double stepSize) throws IllegalStateException, UnsupportedOperationException {
		throw new UnsupportedOperationException();
	}

	@Override
	public void setMinimumStepSize(double minStepSize) throws IllegalStateException, UnsupportedOperationException {
		throw new UnsupportedOperationException();
	}

	@Override
	public void setMaximumStepSize(double maxStepSize) throws IllegalStateException, UnsupportedOperationException {
		throw new UnsupportedOperationException();
	}

}
