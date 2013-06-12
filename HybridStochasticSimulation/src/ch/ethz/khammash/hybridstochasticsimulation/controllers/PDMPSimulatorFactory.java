package ch.ethz.khammash.hybridstochasticsimulation.controllers;

import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPSimulator;

public interface PDMPSimulatorFactory<S extends PDMPSimulator<?>> extends SimulatorFactory<S> {

	public void setIntegratorFactory(IntegratorFactory integratorFactory);

}
