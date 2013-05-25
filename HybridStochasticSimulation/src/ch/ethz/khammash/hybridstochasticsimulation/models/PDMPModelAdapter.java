package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.Collection;
import java.util.Collections;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.events.EventHandler;

public class PDMPModelAdapter implements PDMPModel,
		FirstOrderDifferentialEquations, ReactionNetworkModel, EventHandler {

	private int baseODEDimension;
	private HybridModel hybridModel;
	private double[] propVector;

	public PDMPModelAdapter(PDMPModelAdapter model) {
		this(model.hybridModel);
	}

	public PDMPModelAdapter(FirstOrderDifferentialEquations baseODE, ReactionNetworkModel reactionModel) {
		this(new HybridModelAdapter(baseODE, reactionModel));
	}

	public PDMPModelAdapter(HybridModel hybridModel) {
		this.hybridModel = hybridModel;
		baseODEDimension = hybridModel.getDimension();
		propVector = new double[hybridModel.getPropensityDimension()];
	}

	public FirstOrderDifferentialEquations getBaseODE() {
		return hybridModel;
	}

	@Override
	public boolean hasDeterministicPart() {
		if (hybridModel == null)
			return true;
		else
			return hybridModel.hasDeterministicPart();
	}

	public int getBaseODEDimension() {
		return baseODEDimension;
	}

	@Override
	public int getPropensityDimension() {
		return hybridModel.getPropensityDimension();
	}

	@Override
	public void computePropensities(double t, double[] x, double[] propensities) {
		hybridModel.computePropensities(t, x, propensities);
	}

	@Override
	public void updateState(int reaction, double t, double[] x) {
		hybridModel.updateState(reaction, t, x);
	}

	@Override
	public int getDimension() {
		return baseODEDimension + 2;
	}

	@Override
	public void computeDerivatives(double t, double[] x, double[] xDot) {
		hybridModel.computeDerivatives(t, x, xDot);
		xDot[xDot.length - 2] = 0;
		computePropensities(t, x, propVector);
		for (int i = 0; i < propVector.length; i++)
			xDot[xDot.length - 2] += propVector[i];
		xDot[xDot.length - 1] = 0;
	}

	@Override
	public double g(double t, double[] x) {
		return x[x.length - 1] - x[x.length - 2];
	}

	@Override
	public Action eventOccurred(double t, double[] x, boolean increasing) {
		return Action.STOP;
	}

	@Override
	public void init(double t0, double[] x0, double t) {
	}

	@Override
	public void resetState(double t, double[] x) {
	}

	@Override
	public int getNumberOfSpecies() {
		return getDimension();
	}

	@Override
	public Collection<EventHandler> getOptionalEventHandlers() {
		return Collections.<EventHandler>emptyList();
	}

	@Override
	public FirstOrderDifferentialEquations getFirstOrderDifferentialEquations() {
		return this;
	}

	@Override
	public EventHandler getPDMPEventHandler() {
		return this;
	}

	@Override
	public ReactionNetworkModel getReactionNetworkModel() {
		return this;
	}

	@Override
	public boolean isTimeIndependent() {
		return hybridModel.isTimeIndependent();
	}

	@Override
	public boolean getOptionalEventFlag() {
		return false;
	}

	@Override
	public void handleOptionalEvent(double t, double[] x) {
	}

	@Override
	public void manualCheckOptionalEvent(double t, double[] x) {
	}

	@Override
	public void initialize(double t0, double[] x0) {
	}

}
