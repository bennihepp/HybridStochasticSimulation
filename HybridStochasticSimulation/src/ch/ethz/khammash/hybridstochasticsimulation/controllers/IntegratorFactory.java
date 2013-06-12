package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import org.apache.commons.math3.ode.AbstractIntegrator;

public interface IntegratorFactory {

	public AbstractIntegrator createIntegrator();

}
