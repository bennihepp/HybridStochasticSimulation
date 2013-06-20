package ch.ethz.khammash.hybridstochasticsimulation.factories;

import org.apache.commons.math3.ode.AbstractIntegrator;

public interface IntegratorFactory {

	AbstractIntegrator createIntegrator();

}
