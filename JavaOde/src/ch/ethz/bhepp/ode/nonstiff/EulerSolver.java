package ch.ethz.bhepp.ode.nonstiff;

import ch.ethz.bhepp.ode.EventFunction;
import ch.ethz.bhepp.ode.EventObserver;
import ch.ethz.bhepp.ode.FiniteTimepointProvider;
import ch.ethz.bhepp.ode.FixedStepSolver;
import ch.ethz.bhepp.ode.Ode;
import ch.ethz.bhepp.ode.Solver;
import ch.ethz.bhepp.ode.StateObserver;
import ch.ethz.bhepp.ode.TimepointProvider;
import ch.ethz.bhepp.ode.EventObserver.EventAction;

public class EulerSolver implements Solver, FixedStepSolver {

	public static class EventAtStartOfIntegrationException extends RuntimeException {
		private static final long serialVersionUID = -474930742741047560L;

		public EventAtStartOfIntegrationException(String message) {
			super(message);
		}
	}

    private boolean initialized = false;
	protected Ode ode;
	protected EventFunction ef;
	protected StateObserver stateObserver;
	protected EventObserver eventObserver;
	protected double stepSize;
	protected double[] xDot;
	protected double[] eventValues;
	protected double[] eventValues2;
//	protected int[] eventSigns;
	private TimepointProvider timepointProvider;
	private double[] x;
	private double t;
	private double t1;

	public EulerSolver(double stepSize) {
		this.stepSize = stepSize;
	}

	@Override
	public void initialize(Ode ode) throws InitializationException {
		this.ode = ode;
		this.xDot = new double[ode.getDimensionOfVectorField()];
	}

	@Override
	public void initialize(Ode ode, EventFunction ef,
			StateObserver stateObserver, EventObserver eventObserver) {
		initialize(ode);
		this.ef = ef;
		this.stateObserver = stateObserver;
		this.eventObserver = eventObserver;
		this.eventValues = new double[ef.getNumberOfEventValues()];
		this.eventValues2 = new double[ef.getNumberOfEventValues()];
//		this.eventSigns = new int[ef.getNumberOfEventValues()];
		initialized = true;
	}

	@Override
	public void dispose() {
	}

	@Override
	public void prepareStep(double t0, double[] x0, double t1) throws IllegalStateException {
    	if (ode == null)
    		throw new IllegalStateException("Solver has not yet been initialized");
		this.t = t0;
		this.x = x0;
		this.t1 = t1;
	}

	@Override
	public double integrateStep(double t, double[] x, double t1, double[] xOut) {
		double _stepSize = stepSize;
		if (t + _stepSize > t1)
			_stepSize = t1 - t;
		findNextStateExplicitFixed(t, x, _stepSize, xOut);
		return t + _stepSize;
	}

	@Override
	public double integrateStep(double t, double[] x, double[] xOut, double step) {
		findNextStateExplicitFixed(t, x, step, xOut);
		return t + step;
	}

	@Override
	public double integrateStep() throws IllegalStateException {
		return integrateStep(stepSize);
	}

	@Override
	public double integrateStep(double stepSize) throws IllegalStateException {
		if (t + stepSize > t1)
			stepSize = t1 - t;
		findNextStateExplicitFixed(t, x, stepSize, x);
		t += stepSize;
		return t;
	}

	@Override
	public double integrate(double t0, double[] x0, double t1) {
    	if (!initialized)
    		throw new IllegalStateException("Solver has not yet been initialized");
		double[] x = x0;
		double t = t0;
		ef.computeEventValues(t, x, eventValues);
		for (int i=0; i < ef.getNumberOfEventValues(); i++) {
			if (eventValues[i] == 0)
				throw new EventAtStartOfIntegrationException("Event function has zeros at start of integration");
//			eventSigns[i] = (eventValues[i] > 0) ? 1 : -1;
		}
		stateObserver.report(t, x.clone());
		while (t < t1) {
			double tOld = t;
			double[] xOld = x.clone();
			t = findNextState(t, x, stepSize, t1 - t);
			ef.computeEventValues(t, x, eventValues2);
			for (int k=0; k < ef.getNumberOfEventValues(); k++)
				if (eventValues[k] * eventValues2[k] < 0) {
					// Event occurred during last step, find a better approximation of the time that the event occured
					double tDelta = - (t - tOld) / (eventValues2[k] - eventValues[k]) * eventValues[k];
					x = xOld.clone();
					t = findNextState(tOld, x, tDelta, tDelta);
					if (eventObserver.report(k, t, x) == EventAction.STOP)
						return t;
        		}
			// Switch eventValues and eventValues2 (so eventValues always contains the last evaluation of the event functions)
			double[] tmp = eventValues;
			eventValues = eventValues2;
			eventValues2 = tmp;
			stateObserver.report(t, x.clone());
		}
		return t;
	}

	@Override
	public void prepare(TimepointProvider timepointProvider, double[] x) {
    	if (!initialized)
    		throw new IllegalStateException("Solver has not yet been initialized");
		this.timepointProvider = timepointProvider;
		this.x = x;
	}

	@Override
	public double integrate() {
		double t = timepointProvider.getCurrentTimepoint();
		double tReport = timepointProvider.getNextTimepoint(t);
		ef.computeEventValues(t, x, eventValues);
		for (int i=0; i < ef.getNumberOfEventValues(); i++) {
			if (eventValues[i] == 0)
				throw new EventAtStartOfIntegrationException("Event function has zeros at start of integration");
//			eventSigns[i] = (eventValues[i] > 0) ? 1 : -1;
		}
		while (true) {
			double tOld = t;
			double[] xOld = x.clone();
//			if (Double.isNaN(x[0]))
//				// TODO: Handle NaNs
			t = findNextState(t, x, stepSize, tReport - t);
			ef.computeEventValues(t, x, eventValues2);
			for (int k=0; k < ef.getNumberOfEventValues(); k++)
				if (eventValues[k] * eventValues2[k] < 0) {
					// Event occurred during last step, find a better approximation of the time that the event occured
					double tDelta = - (t - tOld) / (eventValues2[k] - eventValues[k]) * eventValues[k];
					x = xOld.clone();
					t = findNextState(tOld, x, tDelta, tDelta);
					// Just return current state and time
					if (eventObserver.report(k, t, x) == EventAction.STOP)
						return t;
        		}
			// Switch eventValues and eventValues2 (so eventValues always contains the last evaluation of the event functions)
			double[] tmp = eventValues;
			eventValues = eventValues2;
			eventValues2 = tmp;
			if (t == tReport) {
				stateObserver.report(t, x);
	    		if (timepointProvider.hasNextTimepoint(t))
	    			tReport = timepointProvider.getNextTimepoint(t);
	    		else
	    			break;
			}
		}
		return t;
	}

	@Override
	public double integrate(TimepointProvider timepointProvider, double[] x) {
		prepare(timepointProvider, x);
		return integrate();
	}

	protected double findNextState(double t, double[] x, double stepSize, double maxStep) {
		if (stepSize > maxStep)
			stepSize = maxStep;
		findNextStateExplicitFixed(t, x, stepSize, x);
		return t + stepSize;
	}

	protected void findNextStateExplicitFixed(double t, double[] x, double step, double[] out) {
		ode.computeVectorField(t, x, xDot);
		for (int j=0; j < ode.getDimensionOfVectorField(); j++)
			out[j] = x[j] + step * xDot[j];
	}

    public static void main(String args[]) {
    	MyOde ode = new MyOde();
    	MyEf ef = new MyEf();
    	MyObserver observer = new MyObserver();
    	MyEventObserver eventObserver = new MyEventObserver();
    	double stepSize = 0.1;
    	EulerSolver q = new EulerSolver(stepSize);
    	q.initialize(ode, ef, observer, eventObserver);
    	double[] x0 = { 0 };
    	double t0 = 0.0;
    	double t1 = 1.0;
    	q.integrate(t0, x0, t1);
    	x0[0] = 0.0;
    	double[] tSeries = { 0, 0.5, 1.0 };
    	FiniteTimepointProvider timepointProvider = new FiniteTimepointProvider(tSeries);
    	q.integrate(timepointProvider, x0);
        q.dispose();
    };

    protected static class MyOde implements Ode {

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

    protected static class MyEf implements EventFunction {

		@Override
		public int getNumberOfEventValues() {
			return 1;
		}

		@Override
		public void computeEventValues(double t, double[] x, double[] values) {
			values[0] = 2.0 - x[0];
			values[0] = 0.5;
		}
    	
    }

    protected static class MyObserver implements StateObserver {

    	@Override
    	public void report(double t, double[] x) {
    		System.out.println("Report x["+t+"] = "+x[0]);
    	}

		@Override
		public void initialize(double t0, double[] x0, double t1) {
			report(t0, x0);
		}

    }

    protected static class MyEventObserver implements EventObserver {

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
	public void setStepSize(double stepSize) throws UnsupportedOperationException {
		this.stepSize = stepSize;
	}

}
