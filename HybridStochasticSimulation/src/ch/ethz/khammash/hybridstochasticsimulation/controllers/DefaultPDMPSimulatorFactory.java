package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import org.apache.commons.math3.ode.AbstractIntegrator;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPSimulator;

public class DefaultPDMPSimulatorFactory implements PDMPSimulatorFactory {

	private double ehMaxCheckInterval;
	private double ehConvergence;
	private double ehConvergenceFactor;
	private int ehMaxIterationCount;

	public DefaultPDMPSimulatorFactory() {
		ehMaxCheckInterval = PDMPSimulator.DEFAULT_EVENT_HANDLER_MAX_CHECK_INTERVAL;
		ehConvergence = PDMPSimulator.DEFAULT_EVENT_HANDLER_CONVERGENCE;
		ehConvergenceFactor = PDMPSimulator.DEFAULT_EVENT_HANDLER_CONVERGENCE_FACTOR;
		ehMaxIterationCount = PDMPSimulator.DEFAULT_EVENT_HANDLER_MAX_ITERATION_COUNT;
	}

	@Override
	public PDMPSimulator createSimulator(AbstractIntegrator integrator, RandomDataGenerator rdg) {
		PDMPSimulator sim = new PDMPSimulator(integrator, rdg);
		sim.setEventHandlerConvergence(getEventHandlerConvergence());
		sim.setEventHandlerConvergenceFactor(getEventHandlerConvergenceFactor());
		sim.setEventHandlerMaxCheckInterval(getEventHandlerMaxCheckInterval());
		sim.setEventHandlerMaxIterationCount(getEventHandlerMaxIterationCount());
		return sim;
	}

	public double getEventHandlerConvergence() {
		return ehConvergence;
	}

	public void setEventHandlerConvergence(double convergence) {
		this.ehConvergence = convergence;
		if (!Double.isNaN(convergence))
			this.ehConvergenceFactor = Double.NaN;
	}

	public double getEventHandlerConvergenceFactor() {
		return ehConvergenceFactor;
	}

	public void setEventHandlerConvergenceFactor(double convergenceFactor) {
		this.ehConvergenceFactor = convergenceFactor;
		if (!Double.isNaN(convergenceFactor))
			this.ehConvergence = Double.NaN;
	}

	public double getEventHandlerMaxCheckInterval() {
		return ehMaxCheckInterval;
	}

	public void setEventHandlerMaxCheckInterval(double maxCheckInterval) {
		this.ehMaxCheckInterval = maxCheckInterval;
	}

	public int getEventHandlerMaxIterationCount() {
		return ehMaxIterationCount;
	}

	public void setEventHandlerMaxIterationCount(int maxIterationCount) {
		this.ehMaxIterationCount = maxIterationCount;
	}

}