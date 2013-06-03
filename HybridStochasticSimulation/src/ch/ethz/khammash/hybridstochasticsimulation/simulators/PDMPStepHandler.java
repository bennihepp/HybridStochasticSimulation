package ch.ethz.khammash.hybridstochasticsimulation.simulators;

import org.apache.commons.math3.ode.sampling.StepHandler;

import ch.ethz.khammash.hybridstochasticsimulation.models.PDMPModel;


public interface PDMPStepHandler extends StepHandler {

	public void setPDMPModel(PDMPModel model);

	public void reset();

}
