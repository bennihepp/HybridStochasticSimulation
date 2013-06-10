package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.Collection;
import java.util.Collections;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.events.EventHandler;

public class PDMPModelAdapter implements FirstOrderDifferentialEquations, PDMPModel, EventHandler {

	private int baseODEDimension;
	private HybridModel hybridModel;
	private double[] propVector;

	public PDMPModelAdapter(PDMPModelAdapter model) {
		this(model.hybridModel);
	}

	public PDMPModelAdapter(FirstOrderDifferentialEquations baseODE, ReactionNetworkModel reactionModel) {
		if (reactionModel.getStateDimension() != baseODE.getDimension())
			throw new UnsupportedOperationException("Expected reactionModel.getNumberOfStates() == baseODE.getDimension()");
		setModel(baseODE, reactionModel);
	}

	public PDMPModelAdapter(HybridModel hybridModel) {
		setHybridModel(hybridModel);
	}

	public void setModel(FirstOrderDifferentialEquations baseODE, ReactionNetworkModel reactionModel) {
		setHybridModel(new HybridModelAdapter(baseODE, reactionModel));
	}

	public HybridModel getHybridModel() {
		return hybridModel;
	}

	public void setHybridModel(HybridModel hybridModel) {
		if (hybridModel.getStateDimension() != hybridModel.getDimension())
			throw new UnsupportedOperationException("Invalid hybrid model. Expected getNumberOfStates() == getDimension()");
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
	public int getDimension() {
		return baseODEDimension + 2;
	}

	@Override
	public void computeDerivatives(double t, double[] x, double[] xDot) {
		hybridModel.computeDerivatives(t, x, xDot);
		hybridModel.computePropensities(t, x, propVector);
		xDot[xDot.length - 2] = 0.0;
		xDot[xDot.length - 1] = 0.0;
		for (int i = 0; i < propVector.length; i++)
			xDot[xDot.length - 2] += propVector[i];
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
		return hybridModel;
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

	@Override
	public int getStateDimension() {
		return hybridModel.getStateDimension();
	}

	@Override
	public int getPropensityDimension() {
		return hybridModel.getPropensityDimension();
	}

}
