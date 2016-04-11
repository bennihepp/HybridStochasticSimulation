package ch.ethz.bhepp.ode.boost;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;

import ch.ethz.bhepp.ode.BufferOdeAdapter;
import ch.ethz.bhepp.ode.FixedStepSolver;
import ch.ethz.bhepp.ode.Ode;
import ch.ethz.bhepp.ode.Solver.InitializationException;

public class FixedBoostOdeintSolver implements FixedStepSolver {

    private static final int BOOSTSTEPPERTYPE_DORMANDPRINCE5 = 1;
    private static final int BOOSTSTEPPERTYPE_BULIRSCHSTOER = 2;
    private static final int BOOSTSTEPPERTYPE_ROSENBROCK4 = 3;
 
	public enum BoostStepperType {
		DormandPrince5 (BOOSTSTEPPERTYPE_DORMANDPRINCE5),
		BulirschStoer (BOOSTSTEPPERTYPE_BULIRSCHSTOER),
		Rosenbrock4 (BOOSTSTEPPERTYPE_ROSENBROCK4);

		private int type;

		private BoostStepperType(int type) {
			this.type = type;
		}

		public int getType() {
			return type;
		}
	}

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

	private double stepSize;
	private BoostStepperType stepperType;
	private BufferOdeAdapter ode;
	private double relTolerance;
	private double absTolerance;
    private boolean initialized;
	private boolean prepared;
    private long jni_pointer;
	private DoubleBuffer xBuffer;
	private DoubleBuffer xTmpBuffer;
	private DoubleBuffer xDotBuffer;
	private double[] x;

	public FixedBoostOdeintSolver(double stepSize, BoostStepperType type) {
		ensureBoostJNILibraryLoaded();
		initialized = false;
		prepared = false;
		this.stepSize = stepSize;
		this.stepperType = type;
		relTolerance = 1e-3;
		absTolerance = 1e-3;
	}

	private native long jni_initialize(BufferOdeAdapter ode,
			DoubleBuffer xBuffer, DoubleBuffer xTmpBuffer, DoubleBuffer xDotBuffer,
			double stepSize, double relTol, double absTol,
			int stepperType);

    private native void jni_dispose(long jni_pointer);

    private native void jni_prepareStep(long jni_pointer, double t0, double t1);

    private native void jni_reinitializeOneStep(long jni_pointer, double t);

    private native double jni_integrateOneStep(long jni_pointer);

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

    public void setRelativeTolerance(double relTolerance) {
//    	checkArgument(relTol >= 0);
    	this.relTolerance = relTolerance;
    }

    public void setAbsoluteTolerance(double absTolerance) {
//    	checkArgument(absTol >= 0);
    	this.absTolerance = absTolerance;
    }

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
    				stepSize, relTolerance, absTolerance, stepperType.getType());
    	} catch (JniException e) {
    		throw new InitializationException("Failed to initialize native CVodeSolver", e);
    	}
        if (jni_pointer == 0)
        	throw new JniException("Failed to initialize native CVodeSolver");
        initialized = true;
		prepared = false;
	}

	private void dispose() {
    	if (initialized) {
    		xBuffer = null;
    		xTmpBuffer = null;
    		xDotBuffer = null;
    		jni_dispose(jni_pointer);
    		initialized = false;
    	}
	}

	@Override
	public void prepareStep(double t0, double[] x, double t1) throws IllegalStateException {
		this.x = x;
		xBuffer.position(0);
		xBuffer.put(this.x);
		jni_prepareStep(jni_pointer, t0, t1);
	}

	@Override
	public double integrateStep() {
		double tNext = jni_integrateOneStep(jni_pointer);
		xBuffer.position(0);
		xBuffer.get(this.x);
		return tNext;
	}

	@Override
	public double integrateStep(double step) {
		return Double.NaN;
	}

	@Override
	public double integrateStep(double t, double[] x, double t1, double[] xOut) {
		jni_prepareStep(jni_pointer, t, t1);
		return integrateStep();
	}

	@Override
	public double integrateStep(double t, double[] x, double[] xOut, double step) {
		return Double.NaN;
	}

    public static void main(String args[]) {
    	double stepSize = 0.01;
    	FixedBoostOdeintSolver solver = new FixedBoostOdeintSolver(stepSize, FixedBoostOdeintSolver.BoostStepperType.DormandPrince5);
    	MyOde ode = new MyOde();
    	solver.initialize(ode);
    	double t0 = 0.0;
    	double t1 = 1.0;
    	double t = t0;
    	double[] x = new double[ode.getDimensionOfVectorField()];
    	x[0] = 0.2;
    	solver.prepareStep(t0, x, t1);
    	while (t < t1) {
    		t = solver.integrateStep();
    		System.out.println("x(" + t + ")=" + x[0]);
    	}
    	solver.dispose();
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

	@Override
	public void setStepSize(double stepSize) throws UnsupportedOperationException {
		throw new UnsupportedOperationException();
	}

}
