package ch.ethz.khammash.hybridstochasticsimulation.models;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;


public class HybridModelAdapter implements HybridModel {

	private FirstOrderDifferentialEquations deterministicModel;
	private StochasticReactionNetworkModel stochasticModel;

	public HybridModelAdapter(FirstOrderDifferentialEquations deterministicModel, StochasticReactionNetworkModel stochasticModel) {
		if (stochasticModel.getNumberOfSpecies() != deterministicModel.getDimension())
			throw new UnsupportedOperationException("Expected stochasticModel.getNumberOfSpecies() == deterministicModel.getDimension()");
		this.deterministicModel = deterministicModel;
		this.stochasticModel = stochasticModel;
	}

	@Override
	public boolean hasDeterministicPart() {
		return true;
	}

	@Override
	public boolean isTimeIndependent() {
		return false;
	}

	@Override
	public FirstOrderDifferentialEquations getDeterministicModel() {
		return deterministicModel;
	}

	@Override
	public StochasticReactionNetworkModel getStochasticModel() {
		return stochasticModel;
	}

}
