package ch.ethz.khammash.hybridstochasticsimulation.simulators;

import java.util.ArrayList;
import java.util.Collection;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;
import org.apache.commons.math3.ode.sampling.StepHandler;

import ch.ethz.khammash.hybridstochasticsimulation.models.DeterministicModel;


public class DeterministicModelSimulator {

	protected FirstOrderIntegrator integrator;
	protected Collection<StepHandler> stepHandlers;

	public DeterministicModelSimulator() {
		this(null);
	}

	public DeterministicModelSimulator(FirstOrderIntegrator integrator) {
		if (integrator == null)
			integrator = new DormandPrince853Integrator(1.0e-8, 100.0, 1.0e-10, 1.0e-10);
		else
			this.integrator = integrator;
		stepHandlers = new ArrayList<StepHandler>();
	}

	public double simulate(DeterministicModel model, final double t0, final double[] x0, double t1, double[] x1) {
    	FirstOrderDifferentialEquations ode = model.getFirstOrderDifferentialEquations();
		double[] x = new double[x0.length];
		for (int i=0; i < x0.length; i++)
			x[i] = x0[i];
		for (StepHandler handler : stepHandlers) {
			handler.init(t0, x, t1);
			integrator.addStepHandler(handler);
		}
//    	long evaluationCounter = 0;
		// Evolve ODE until end timepoint
    	double t = integrator.integrate(ode, t0, x, t1, x);
//    	evaluationCounter += integrator.getEvaluations();
//		System.out.println("Total of " + evaluationCounter + " evaluations");
		integrator.clearStepHandlers();
		for (int i=0; i < x1.length; i++)
			x1[i] = x[i];
		return t;
	}

	public void addStepHandler(StepHandler handler) {
		if (stepHandlers.contains(handler) == false)
			stepHandlers.add(handler);
	}
	public void removeStepHandler(StepHandler handler) {
		if (stepHandlers.contains(handler))
			stepHandlers.remove(handler);
	}
	public Collection<StepHandler> getStepHandlers() {
		return stepHandlers;
	}
	public void clearStepHandlers() {
		stepHandlers.clear();
	}

}
