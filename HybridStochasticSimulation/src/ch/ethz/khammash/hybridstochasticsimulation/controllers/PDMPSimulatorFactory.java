package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import org.apache.commons.math3.ode.AbstractIntegrator;
import org.apache.commons.math3.random.RandomDataGenerator;

import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPSimulator;

public interface PDMPSimulatorFactory {
	public PDMPSimulator createSimulator(AbstractIntegrator integrator, RandomDataGenerator rdg);
}