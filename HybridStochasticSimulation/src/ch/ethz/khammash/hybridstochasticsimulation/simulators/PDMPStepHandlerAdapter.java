package ch.ethz.khammash.hybridstochasticsimulation.simulators;

import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;

public class PDMPStepHandlerAdapter implements PDMPStepHandler {

	StepHandler handler;
	boolean initialized;

	public PDMPStepHandlerAdapter(StepHandler handler) {
		this.handler = handler;
		initialized = false;
	}

	@Override
	public void reset() {
		initialized = false;
	}

	@Override
	public void setPDMPModel(PDMPModel model) {
	}

	@Override
	public void handleStep(StepInterpolator interpolator, boolean isLast) throws MaxCountExceededException {
		handler.handleStep(interpolator, isLast);
	}

	@Override
	public void init(double t0, double[] y0, double t) {
		if (initialized == false) {
			initialized = true;
			handler.init(t0, y0, t);
		}
	}

}