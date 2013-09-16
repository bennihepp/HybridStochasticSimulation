package ch.ethz.khammash.ode.nonstiff;

import ch.ethz.khammash.ode.EventFunction;
import ch.ethz.khammash.ode.EventObserver;
import ch.ethz.khammash.ode.EventObserver.EventAction;
import ch.ethz.khammash.ode.FiniteTimepointProvider;
import ch.ethz.khammash.ode.Ode;
import ch.ethz.khammash.ode.Solver;
import ch.ethz.khammash.ode.StateObserver;
import ch.ethz.khammash.ode.TimepointProvider;

public class EulerSolver implements Solver {

	public static class EventAtStartOfIntegrationException extends RuntimeException {
		private static final long serialVersionUID = -474930742741047560L;

		public EventAtStartOfIntegrationException(String message) {
			super(message);
		}
	}

	protected Ode ode;
	protected EventFunction ef;
	protected StateObserver stateObserver;
	protected EventObserver eventObserver;
	protected double step;
	protected double[] xDot;
	protected double[] eventValues;
	protected double[] eventValues2;
//	protected int[] eventSigns;
	private TimepointProvider timepointProvider;
	private double[] x;

	public EulerSolver(double step) {
		this.step = step;
	}

	@Override
	public void initialize(Ode ode, EventFunction ef,
			StateObserver stateObserver, EventObserver eventObserver) {
		this.ode = ode;
		this.ef = ef;
		this.stateObserver = stateObserver;
		this.eventObserver = eventObserver;
		this.xDot = new double[ode.getDimensionOfVectorField()];
		this.eventValues = new double[ef.getNumberOfEventValues()];
		this.eventValues2 = new double[ef.getNumberOfEventValues()];
//		this.eventSigns = new int[ef.getNumberOfEventValues()];
	}

	@Override
	public void dispose() {
	}

	@Override
	public double integrate(double t0, double[] x0, double t1) {
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
			t = findNextState(t, x, step, t1 - t);
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
			t = findNextState(t, x, step, tReport - t);
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

	protected double findNextState(double t, double[] x, double step, double maxStep) {
		if (step > maxStep)
			step = maxStep;
		findNextStateExplicitFixed(t, x, step, x);
		return t + step;
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
    	double step = 0.1;
    	EulerSolver q = new EulerSolver(step);
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

}
