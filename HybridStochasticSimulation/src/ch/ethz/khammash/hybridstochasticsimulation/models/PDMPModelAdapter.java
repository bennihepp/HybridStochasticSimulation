package ch.ethz.khammash.hybridstochasticsimulation.models;

import java.util.Collections;
import java.util.List;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import ch.ethz.khammash.hybridstochasticsimulation.simulators.PDMPEventObserver;

public class PDMPModelAdapter<T extends HybridModel> implements PDMPModel, FirstOrderDifferentialEquations, PDMPEventObserver {

	private T hybridModel;
	private StochasticReactionNetworkModel transitionMeasure;
	private FirstOrderDifferentialEquations vectorField;
	private double[] propVector;

	public PDMPModelAdapter(PDMPModelAdapter<T> model) {
		this(model.hybridModel);
	}

	public PDMPModelAdapter(T hybridModel) {
		setHybridModel(hybridModel);
	}

	public T getHybridModel() {
		return hybridModel;
	}

	public void setHybridModel(T hybridModel) {
		transitionMeasure = hybridModel.getStochasticModel();
		vectorField = hybridModel.getDeterministicModel();
		this.hybridModel = hybridModel;
		propVector = new double[transitionMeasure.getNumberOfReactions()];
	}

	@Override
	public int getNumberOfSpecies() {
		return hybridModel.getNumberOfSpecies();
	}

	@Override
	public int getNumberOfReactions() {
		return hybridModel.getNumberOfReactions();
	}

	@Override
	public boolean hasVectorField() {
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
	public FirstOrderDifferentialEquations getVectorField() {
		return this;
	}

	@Override
	public StochasticReactionNetworkModel getTransitionMeasure() {
		return transitionMeasure;
	}

	@Override
	public List<PDMPEventObserver> getOptionalEventObservers() {
		return Collections.<PDMPEventObserver>emptyList();
	}

	@Override
	public PDMPEventObserver getJumpEventObserver() {
		return this;
	}

	@Override
	public boolean hasOptionalEventOccured() {
		return false;
	}

	@Override
	public void handleOptionalEvent(double t, double[] x) {
	}

	@Override
	public void checkAndHandleOptionalEvent(double t, double[] x) {
	}

	@Override
	public void initialize(double t0, double[] x0) {
	}

	@Override
	public int getDimension() {
		return vectorField.getDimension() + 2;
	}

	@Override
	public void computeDerivatives(double t, double[] x, double[] xDot) {
		vectorField.computeDerivatives(t, x, xDot);
		transitionMeasure.computePropensities(t, x, propVector);
		xDot[xDot.length - 2] = 0.0;
		xDot[xDot.length - 1] = 0.0;
		for (int i = 0; i < propVector.length; i++)
			xDot[xDot.length - 2] += propVector[i];
	}

	@Override
	public void init(double t0, double[] x0, double t) {
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
	public void resetState(double t, double[] x) {
	}

}
