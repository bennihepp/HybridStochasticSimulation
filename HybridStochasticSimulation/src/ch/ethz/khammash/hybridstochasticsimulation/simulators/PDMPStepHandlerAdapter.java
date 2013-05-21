package ch.ethz.khammash.hybridstochasticsimulation.simulators;

import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

public class PDMPStepHandlerAdapter implements StepHandler {
	StepHandler handler;
	boolean initialized;

	public PDMPStepHandlerAdapter(StepHandler handler) {
		this.handler = handler;
		initialized = false;
	}

	public void reset() {
		initialized = false;
	}

	public void handleStep(StepInterpolator interpolator, boolean isLast) throws MaxCountExceededException {
		handler.handleStep(interpolator, isLast);
	}

	public void init(double t0, double[] y0, double t) {
		if (initialized == false) {
			initialized = true;
			handler.init(t0, y0, t);
		}
	}
}