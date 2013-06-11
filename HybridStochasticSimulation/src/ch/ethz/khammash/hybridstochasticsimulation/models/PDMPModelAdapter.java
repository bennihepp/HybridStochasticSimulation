package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.Collection;
import java.util.Collections;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.events.EventHandler;

public class PDMPModelAdapter implements PDMPModel, FirstOrderDifferentialEquations, EventHandler {

	private HybridModel hybridModel;
	private StochasticReactionNetworkModel stochasticModel;
	private FirstOrderDifferentialEquations deterministicModel;
	private double[] propVector;

	public PDMPModelAdapter(PDMPModelAdapter model) {
		this(model.hybridModel);
	}

	public PDMPModelAdapter(HybridModel hybridModel) {
		setHybridModel(hybridModel);
	}

	public HybridModel getHybridModel() {
		return hybridModel;
	}

	public void setHybridModel(HybridModel hybridModel) {
		stochasticModel = hybridModel.getStochasticModel();
		deterministicModel = hybridModel.getDeterministicModel();
		this.hybridModel = hybridModel;
		propVector = new double[stochasticModel.getNumberOfReactions()];
	}

	@Override
	public boolean hasDeterministicPart() {
		if (hybridModel == null)
			return true;
		else
			return hybridModel.hasDeterministicPart();
	}

	@Override
	public boolean isTimeIndependent() {
		return hybridModel.isTimeIndependent();
	}

	@Override
	public FirstOrderDifferentialEquations getDeterministicModel() {
		return this;
	}

	@Override
	public StochasticReactionNetworkModel getStochasticModel() {
		return stochasticModel;
	}

	@Override
	public Collection<EventHandler> getOptionalEventHandlers() {
		return Collections.<EventHandler>emptyList();
	}

	@Override
	public EventHandler getPDMPEventHandler() {
		return this;
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
	public int getDimension() {
		return deterministicModel.getDimension() + 2;
	}

	@Override
	public void computeDerivatives(double t, double[] x, double[] xDot) {
		deterministicModel.computeDerivatives(t, x, xDot);
		stochasticModel.computePropensities(t, x, propVector);
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

}
