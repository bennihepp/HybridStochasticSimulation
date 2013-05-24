package ch.ethz.khammash.hybridstochasticsimulation.models;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

public class HybridModelAdapter implements HybridModel {

	private FirstOrderDifferentialEquations baseODE;
	private ReactionNetworkModel reactionModel;

	public HybridModelAdapter(FirstOrderDifferentialEquations baseODE, ReactionNetworkModel reactionModel) {
		this.baseODE = baseODE;
		this.reactionModel = reactionModel;
	}

	@Override
	public void computeDerivatives(double t, double[] x, double[] xDot)
			throws MaxCountExceededException, DimensionMismatchException {
		baseODE.computeDerivatives(t, x, xDot);
	}

	@Override
	public int getDimension() {
		return baseODE.getDimension();
	}

	@Override
	public int getNumberOfSpecies() {
		return reactionModel.getNumberOfSpecies();
	}

	@Override
	public int getPropensityDimension() {
		return reactionModel.getPropensityDimension();
	}

	@Override
	public void computePropensities(double t, double[] x, double[] propensities) {
		reactionModel.computePropensities(t, x, propensities);
	}

	@Override
	public void updateState(int reaction, double t, double[] x) {
		reactionModel.updateState(reaction, t, x);
	}

	@Override
	public boolean hasDeterministicPart() {
		return true;
	}

	@Override
	public boolean isTimeIndependent() {
		return false;
	}

}
